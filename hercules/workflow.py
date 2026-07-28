"""Executable HERCULES orchestration around the recovered scientific runtime."""

from __future__ import annotations

import json
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .config import HerculesConfig
from .exceptions import HerculesError
from .m3 import M3_SCIENTIFIC_MODEL, integrate_posterior_tables
from .manifest import RunManifest, checkpoint_matches, checkpoint_path, configuration_hash
from .process import run_process
from .resources import ensemble_script_path
from .scoring import aggregate_score_files, run_plink2_score
from .stages import STAGES, StageSpec, execution_order
from .sumstats import prepare_fastgwa_sumstats


STAGE3_SCIENTIFIC_MODEL = "two-score-validation-trained-superlearner-v1"


@dataclass(frozen=True, slots=True)
class WorkflowPlan:
    trait: str
    configuration_hash: str
    stages: tuple[StageSpec, ...]
    output_dir: Path

    def as_dict(self) -> dict[str, Any]:
        return {
            "trait": self.trait,
            "configuration_hash": self.configuration_hash,
            "stages": [stage.stage_id for stage in self.stages],
            "output_dir": str(self.output_dir),
        }


def plan_workflow(config: HerculesConfig, target: str = "ensemble") -> WorkflowPlan:
    stages = execution_order(target)
    return WorkflowPlan(
        trait=config.trait_name,
        configuration_hash=_configuration_hash(config),
        stages=stages,
        output_dir=Path(config.output_dir),
    )


def execute_workflow(config: HerculesConfig, target: str = "ensemble") -> None:
    """Execute M1/M2/M3 and ensemble stages through *target*.

    M1 and M2 invoke the installed HERCULES inference module. M3 retains the
    validated integration equations, including the current variance-column
    interpretation. R is used for the final SuperLearner procedure.

    The workflow is designed to preserve identified scientific defaults, but
    scientific equivalence is only established for components explicitly
    reported in ``docs/BASELINE_BEHAVIOR.md``.
    """

    output_dir = Path(config.output_dir).resolve()
    temporary_dir = Path(config.temporary_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_dir.mkdir(parents=True, exist_ok=True)
    manifest = create_manifest(config, target)
    manifest.write(output_dir / "run-manifest.json")

    for stage in execution_order(target):
        if stage.stage_id in {"m1", "m2"}:
            _run_inference_stage(config, stage.stage_id)
        elif stage.stage_id == "m3":
            _run_m3(config)
        elif stage.stage_id == "ensemble":
            _run_ensemble(config)
        else:  # pragma: no cover - guarded by the stage registry
            raise HerculesError(f"Unsupported workflow stage: {stage.stage_id}")


def create_manifest(config: HerculesConfig, target: str = "ensemble") -> RunManifest:
    plan = plan_workflow(config, target)
    return RunManifest.create(
        version=__version__,
        trait=config.trait_name,
        config=config.as_dict(),
        stages=tuple(stage.stage_id for stage in plan.stages),
    )


def _run_inference_stage(config: HerculesConfig, stage_id: str) -> None:
    stage = STAGES[stage_id]
    model_config = config.m1 if stage_id == "m1" else config.m2
    ancestry = config.target_ancestry if stage_id == "m1" else config.base_ancestry
    sumstats = (
        config.inputs.summary_statistics.target_path
        if stage_id == "m1"
        else config.inputs.summary_statistics.base_path
    )
    ld_reference = config.inputs.ld_reference["target" if stage_id == "m1" else "base"]
    validation_genotype = str(
        model_config.get(
            "validation_genotype",
            config.inputs.validation_genotype
            if stage_id == "m1"
            else config.inputs.genotype_prefixes.get(
                "base_validation", config.inputs.validation_genotype
            ),
        )
    )
    validation_phenotype = str(model_config.get("validation_phenotype", config.inputs.phenotype_file))
    validation_keep = str(model_config.get("validation_keep", ""))
    config_hash = _configuration_hash(config)
    selected_posteriors: list[Path] = []
    hyperparameter_tables: list[Path] = []
    selection_tables: list[Path] = []
    for chromosome in config.chromosomes:
        marker = checkpoint_path(
            config.output_dir,
            stage=stage,
            trait=config.trait_name,
            chromosome=chromosome,
        )
        stage_prefix = Path(config.temporary_dir) / f"{stage.output_prefix}.{chromosome}"
        selected_posterior = Path(f"{stage_prefix}.fit.gz")
        hyperparameters = Path(f"{stage_prefix}.hyp")
        selection = Path(f"{stage_prefix}.validation")

        complete = (
            checkpoint_matches(marker, config_hash)
            and selected_posterior.is_file()
            and hyperparameters.is_file()
            and selection.is_file()
        )
        if not (config.checkpoint.enabled and config.checkpoint.resume and complete):
            sumstats_format = str(model_config.get("sumstats_format", "fastgwa")).lower()
            source_sumstats = _format_path(
                sumstats, chromosome=chromosome, ancestry=ancestry
            )
            if sumstats_format == "fastgwa":
                prepared_sumstats = prepare_fastgwa_sumstats(
                    source_sumstats,
                    Path(config.temporary_dir)
                    / f"{stage.output_prefix}.{chromosome}.internal.fastGWA.tsv",
                )
                inference_sumstats = prepared_sumstats.path
            else:
                inference_sumstats = Path(source_sumstats)
            command: list[str | Path] = [
                *_fit_command_prefix(),
                "-l",
                _format_path(ld_reference, chromosome=chromosome, ancestry=ancestry),
                "-s",
                inference_sumstats,
                "--output-dir",
                config.temporary_dir,
                "--output-file-prefix",
                f"{stage.output_prefix}.{chromosome}",
                "--hyp-search",
                _hyperparameter_search_code(model_config),
                "--pi-steps",
                str(model_config.get("pi_steps", 10)),
                "--sigma-epsilon-steps",
                str(model_config.get("sigma_epsilon_steps", 10)),
                "--sumstats-format",
                sumstats_format,
                "--backend",
                str(model_config.get("backend", "plink")),
                "--validation-bfile",
                _format_path(validation_genotype, chromosome=chromosome, ancestry=ancestry),
                "--validation-pheno",
                _format_path(validation_phenotype, chromosome=chromosome, ancestry=ancestry),
                "--temp-dir",
                config.temporary_dir,
                "--grid-metric",
                str(model_config.get("grid_metric", "validation")),
                "--validation-metric",
                "auc" if config.inputs.trait_type == "binary" else "r2",
                "--max-iter",
                str(model_config.get("max_iter", 500)),
                "--threads",
                str(config.execution.threads),
                "--n-jobs",
                str(config.execution.parallel_jobs),
                "--seed",
                str(config.execution.seed),
            ]
            if validation_keep:
                command.extend(
                    [
                        "--validation-keep",
                        _format_path(validation_keep, chromosome=chromosome, ancestry=ancestry),
                    ]
                )
            command.extend(_per_snp_prior_cli_args(model_config))
            result = run_process(command, env=_scientific_subprocess_environment(config))
            Path(f"{stage_prefix}.runner.stdout.log").write_text(result.stdout, encoding="utf-8")
            Path(f"{stage_prefix}.runner.stderr.log").write_text(result.stderr, encoding="utf-8")
            if not selected_posterior.is_file() or not hyperparameters.is_file() or not selection.is_file():
                raise HerculesError(
                    f"{stage.display_name} did not produce expected posterior files for chromosome "
                    f"{chromosome}: {selected_posterior}, {hyperparameters}, {selection}"
                )
            _write_checkpoint(marker, config_hash, stage_id, chromosome)

        selected_posteriors.append(selected_posterior)
        hyperparameter_tables.append(hyperparameters)
        selection_tables.append(selection)

    _concatenate_tables(
        selected_posteriors,
        Path(config.output_dir) / f"{stage.output_prefix}.selected-posterior.tsv.gz",
    )
    _concatenate_tables(
        hyperparameter_tables,
        Path(config.output_dir) / f"{stage.output_prefix}.selected-hyperparameters.tsv",
    )
    _concatenate_tables(
        selection_tables,
        Path(config.output_dir) / f"{stage.output_prefix}.selection-metrics.tsv",
    )


def _run_m3(config: HerculesConfig) -> None:
    stage = STAGES["m3"]
    config_hash = _configuration_hash(config)
    posterior_paths: list[Path] = []
    diagnostic_paths: list[Path] = []
    for chromosome in config.chromosomes:
        marker = checkpoint_path(
            config.output_dir,
            stage=stage,
            trait=config.trait_name,
            chromosome=chromosome,
        )
        m1 = Path(config.temporary_dir) / f"{STAGES['m1'].output_prefix}.{chromosome}.fit.gz"
        m2 = Path(config.temporary_dir) / f"{STAGES['m2'].output_prefix}.{chromosome}.fit.gz"
        m3 = Path(config.temporary_dir) / f"{stage.output_prefix}.{chromosome}.posterior.tsv"
        diagnostics = Path(config.temporary_dir) / f"{stage.output_prefix}.{chromosome}.diagnostics.tsv"
        complete = (
            checkpoint_matches(marker, config_hash)
            and m3.is_file()
            and diagnostics.is_file()
        )
        if not (config.checkpoint.enabled and config.checkpoint.resume and complete):
            if not m1.is_file() or not m2.is_file():
                raise HerculesError(f"M3 requires M1/M2 posterior files: {m1}, {m2}")
            integrate_posterior_tables(
                m1,
                m2,
                m3,
                diagnostics_path=diagnostics,
                max_iter=int(config.m3.get("max_iter", 1000)),
                tol=float(config.m3.get("tol", 1e-6)),
                quadrature_points=int(config.m3.get("quadrature_points", 32)),
            )
            _write_checkpoint(marker, config_hash, "m3", chromosome)
        posterior_paths.append(m3)
        diagnostic_paths.append(diagnostics)

    _concatenate_tables(
        posterior_paths,
        Path(config.output_dir) / f"{stage.output_prefix}.calibrated-posterior.tsv.gz",
    )
    _concatenate_tables(
        diagnostic_paths,
        Path(config.output_dir) / f"{stage.output_prefix}.convergence-diagnostics.tsv",
    )


def _run_ensemble(config: HerculesConfig) -> None:
    stage = STAGES["ensemble"]
    config_hash = _configuration_hash(config)
    marker = checkpoint_path(
        config.output_dir,
        stage=stage,
        trait=config.trait_name,
        chromosome=None,
    )
    output_prefix = Path(config.output_dir) / stage.output_prefix
    predictions = Path(f"{output_prefix}.predictions.tsv")
    metrics = Path(f"{output_prefix}.metrics.tsv")
    model_path = Path(f"{output_prefix}.model.rds")
    if (
        config.checkpoint.enabled
        and config.checkpoint.resume
        and checkpoint_matches(marker, config_hash)
        and predictions.is_file()
        and model_path.is_file()
        and (not config.inputs.target_test_phenotype or metrics.is_file())
    ):
        return

    validation_predictors = _build_stage3_score_matrix(
        config,
        cohort="validation",
        genotype=config.inputs.target_validation_genotype,
    )
    test_predictors = _build_stage3_score_matrix(
        config,
        cohort="test",
        genotype=config.inputs.target_test_genotype,
    )
    _ensure_disjoint_iids(validation_predictors, test_predictors)
    validation_path = Path(config.output_dir) / f"{stage.output_prefix}.validation-scores.tsv"
    test_path = Path(config.output_dir) / f"{stage.output_prefix}.test-scores.tsv"
    validation_predictors.to_csv(validation_path, sep="\t", index=False)
    test_predictors.to_csv(test_path, sep="\t", index=False)

    result = run_process(
        [
            config.tools.rscript,
            ensemble_script_path(),
            validation_path,
            config.inputs.target_validation_phenotype,
            test_path,
            config.inputs.target_test_phenotype,
            config.inputs.phenotype_column,
            ",".join(config.inputs.covariates),
            config.inputs.trait_type,
            output_prefix,
            str(config.execution.seed),
        ]
    )
    Path(f"{output_prefix}.stdout.log").write_text(result.stdout, encoding="utf-8")
    Path(f"{output_prefix}.stderr.log").write_text(result.stderr, encoding="utf-8")
    if not predictions.is_file() or not model_path.is_file():
        raise HerculesError("Stage 3 completed without producing predictions and a fitted model")
    if config.inputs.target_test_phenotype and not metrics.is_file():
        raise HerculesError("Stage 3 completed without producing the requested test metric")
    _write_checkpoint(marker, config_hash, "ensemble", None)


def _build_stage3_score_matrix(
    config: HerculesConfig,
    *,
    cohort: str,
    genotype: str,
) -> pd.DataFrame:
    """Score the two manuscript predictors in one target-population cohort."""

    target_scores: list[Path] = []
    calibrated_scores: list[Path] = []
    for chromosome in config.chromosomes:
        genotype_prefix = _format_path(
            genotype,
            chromosome=chromosome,
            ancestry=config.target_ancestry,
        )
        target_posterior = (
            Path(config.temporary_dir)
            / f"{STAGES['m1'].output_prefix}.{chromosome}.fit.gz"
        )
        calibrated_posterior = (
            Path(config.temporary_dir)
            / f"{STAGES['m3'].output_prefix}.{chromosome}.posterior.tsv"
        )
        for label, posterior, collector in (
            ("target_stage1", target_posterior, target_scores),
            ("calibrated_stage2", calibrated_posterior, calibrated_scores),
        ):
            if not posterior.is_file():
                raise HerculesError(f"Stage 3 scoring posterior is missing: {posterior}")
            collector.append(
                run_plink2_score(
                    plink2=config.tools.plink2,
                    genotype_prefix=genotype_prefix,
                    score_file=posterior,
                    output_prefix=(
                        Path(config.temporary_dir)
                        / f"HERCULES_stage3.{cohort}.{label}.{chromosome}"
                    ),
                    grid=False,
                )
            )

    target = _aggregate_single_score(
        target_scores,
        Path(config.temporary_dir) / f"HERCULES_stage3.{cohort}.target.tsv",
        "target_stage1_score",
    )
    calibrated = _aggregate_single_score(
        calibrated_scores,
        Path(config.temporary_dir) / f"HERCULES_stage3.{cohort}.calibrated.tsv",
        "calibrated_stage2_score",
    )
    target_iids = set(target["IID"])
    calibrated_iids = set(calibrated["IID"])
    if target_iids != calibrated_iids:
        target_only = sorted(str(value) for value in target_iids - calibrated_iids)
        calibrated_only = sorted(str(value) for value in calibrated_iids - target_iids)
        raise HerculesError(
            f"Stage 3 {cohort} target/calibrated score sample sets differ; "
            f"target_only={len(target_only)}"
            f"{f' (example {target_only[0]})' if target_only else ''}, "
            f"calibrated_only={len(calibrated_only)}"
            f"{f' (example {calibrated_only[0]})' if calibrated_only else ''}"
        )
    combined = target.merge(
        calibrated,
        on="IID",
        how="inner",
        sort=False,
        validate="one_to_one",
    )
    if combined.empty:
        raise HerculesError(f"Stage 3 {cohort} score intersection is empty")
    return combined.loc[
        :, ["IID", "target_stage1_score", "calibrated_stage2_score"]
    ]


def _aggregate_single_score(
    paths: list[Path],
    intermediate_path: Path,
    output_column: str,
) -> pd.DataFrame:
    aggregate_score_files(paths, intermediate_path)
    table = pd.read_csv(intermediate_path, sep="\t")
    score_columns = [column for column in table.columns if column != "IID"]
    if len(score_columns) != 1:
        raise HerculesError(
            f"Expected exactly one score column for {output_column}; found {score_columns}"
        )
    return table.rename(columns={score_columns[0]: output_column})


def _ensure_disjoint_iids(
    validation_predictors: pd.DataFrame,
    test_predictors: pd.DataFrame,
) -> None:
    overlap = set(validation_predictors["IID"]).intersection(test_predictors["IID"])
    if overlap:
        example = sorted(str(value) for value in overlap)[0]
        raise HerculesError(
            "Target validation and test samples must be disjoint; "
            f"overlapping IID detected: {example}"
        )


def _concatenate_tables(paths: list[Path], output_path: Path) -> Path:
    if not paths:
        raise HerculesError(f"No tables were provided for {output_path}")
    tables = [pd.read_csv(path, sep="\t") for path in paths]
    expected = list(tables[0].columns)
    for path, table in zip(paths[1:], tables[1:]):
        if list(table.columns) != expected:
            raise HerculesError(f"Table schemas differ while creating {output_path}: {path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(output_path, sep="\t", index=False)
    return output_path


def _configuration_hash(config: HerculesConfig) -> str:
    return configuration_hash(
        {
            "package_version": __version__,
            "stage2_scientific_model": M3_SCIENTIFIC_MODEL,
            "stage3_scientific_model": STAGE3_SCIENTIFIC_MODEL,
            "configuration": config.as_dict(),
        }
    )


def _format_path(value: str, *, chromosome: str, ancestry: str) -> str:
    try:
        return value.format(chrom=chromosome, chromosome=chromosome, ancestry=ancestry)
    except KeyError as exc:
        raise HerculesError(f"Unknown placeholder in configured path {value!r}: {exc}") from exc


def _hyperparameter_search_code(config: dict[str, Any]) -> str:
    value = str(config.get("hyperparameter_search", "grid")).lower()
    return {"grid": "GS", "gs": "GS", "em": "EM", "bma": "BMA", "bo": "BO"}.get(
        value, value.upper()
    )


def _per_snp_prior_cli_args(model_config: dict[str, Any]) -> list[str]:
    """Translate the stage-specific tau_beta initialization for fit_cli."""

    prior = model_config.get("per_snp_prior", {})
    if not prior.get("enabled", True):
        return []
    return [
        "--initial-per-snp-prior-column",
        str(prior.get("column", "PVAL")),
        "--per-snp-prior-input-type",
        str(prior.get("input_type", "variance")),
    ]


def _fit_command_prefix() -> tuple[str, ...]:
    """Run the installed module without allowing the checkout to shadow it."""

    return (sys.executable, "-P", "-m", "hercules.fit_cli")


def _scientific_subprocess_environment(config: HerculesConfig) -> dict[str, str]:
    """Expose configured genotype tools to nested magenpy subprocesses."""

    environment = os.environ.copy()
    configured_directories: list[str] = []
    for executable in (config.tools.plink, config.tools.plink2):
        resolved = shutil.which(executable) or executable
        path = Path(resolved).expanduser()
        if path.is_file():
            directory = str(path.resolve().parent)
            if directory not in configured_directories:
                configured_directories.append(directory)
    existing_path = environment.get("PATH", "")
    environment["PATH"] = os.pathsep.join(
        [*configured_directories, *([existing_path] if existing_path else [])]
    )
    return environment


def _write_checkpoint(path: Path, config_hash: str, stage_id: str, chromosome: str | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "status": "complete",
                "configuration_hash": config_hash,
                "stage": stage_id,
                "chromosome": chromosome,
                "package_version": __version__,
                "stage2_scientific_model": M3_SCIENTIFIC_MODEL,
                "stage3_scientific_model": STAGE3_SCIENTIFIC_MODEL,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
