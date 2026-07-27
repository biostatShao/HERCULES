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
from .m3 import integrate_posterior_tables
from .manifest import RunManifest, checkpoint_matches, checkpoint_path, configuration_hash
from .process import run_process
from .resources import ensemble_script_path
from .scoring import aggregate_score_files, run_plink2_score
from .stages import STAGES, StageSpec, execution_order


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
        configuration_hash=configuration_hash(config.as_dict()),
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
    target_genotype = config.inputs.genotype_prefixes.get(
        "target", config.inputs.validation_genotype
    )
    if not target_genotype:
        raise HerculesError("inputs.genotype_prefixes.target is required for PRS scoring")

    config_hash = configuration_hash(config.as_dict())
    selected_scores: list[Path] = []
    grid_scores: list[Path] = []
    for chromosome in config.chromosomes:
        marker = checkpoint_path(
            config.output_dir,
            stage=stage,
            trait=config.trait_name,
            chromosome=chromosome,
        )
        stage_prefix = Path(config.temporary_dir) / f"{stage.output_prefix}.{chromosome}"
        posterior_grid = Path(f"{stage_prefix}.beta.gz")
        selected_posterior = Path(f"{stage_prefix}.fit.gz")
        selected_score_prefix = Path(config.temporary_dir) / f"{stage.output_prefix}.selected.{chromosome}"
        grid_score_prefix = Path(config.temporary_dir) / f"{stage.output_prefix}.grid.{chromosome}"
        selected_score = Path(f"{selected_score_prefix}.sscore")
        grid_score = Path(f"{grid_score_prefix}.sscore")

        complete = (
            checkpoint_matches(marker, config_hash)
            and posterior_grid.is_file()
            and selected_posterior.is_file()
            and selected_score.is_file()
            and grid_score.is_file()
        )
        if not (config.checkpoint.enabled and config.checkpoint.resume and complete):
            command: list[str | Path] = [
                *_fit_command_prefix(),
                "-l",
                _format_path(ld_reference, chromosome=chromosome, ancestry=ancestry),
                "-s",
                _format_path(sumstats, chromosome=chromosome, ancestry=ancestry),
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
                str(model_config.get("sumstats_format", "fastgwa")),
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
            result = run_process(command, env=_scientific_subprocess_environment(config))
            Path(f"{stage_prefix}.runner.stdout.log").write_text(result.stdout, encoding="utf-8")
            Path(f"{stage_prefix}.runner.stderr.log").write_text(result.stderr, encoding="utf-8")
            if not posterior_grid.is_file() or not selected_posterior.is_file():
                raise HerculesError(
                    f"{stage.display_name} did not produce expected posterior files for chromosome "
                    f"{chromosome}: {posterior_grid}, {selected_posterior}"
                )

            genotype = _format_path(target_genotype, chromosome=chromosome, ancestry=config.target_ancestry)
            grid_score = run_plink2_score(
                plink2=config.tools.plink2,
                genotype_prefix=genotype,
                score_file=posterior_grid,
                output_prefix=grid_score_prefix,
                grid=True,
            )
            selected_score = run_plink2_score(
                plink2=config.tools.plink2,
                genotype_prefix=genotype,
                score_file=selected_posterior,
                output_prefix=selected_score_prefix,
                grid=False,
            )
            _write_checkpoint(marker, config_hash, stage_id, chromosome)

        selected_scores.append(selected_score)
        grid_scores.append(grid_score)

    aggregate_score_files(
        selected_scores,
        Path(config.output_dir) / f"{stage.output_prefix}.scores.tsv",
    )
    aggregate_score_files(
        grid_scores,
        Path(config.output_dir) / f"{stage.output_prefix}.grid-scores.tsv",
    )


def _run_m3(config: HerculesConfig) -> None:
    stage = STAGES["m3"]
    config_hash = configuration_hash(config.as_dict())
    grid_scores: list[Path] = []
    target_genotype = config.inputs.genotype_prefixes.get(
        "target", config.inputs.validation_genotype
    )
    for chromosome in config.chromosomes:
        marker = checkpoint_path(
            config.output_dir,
            stage=stage,
            trait=config.trait_name,
            chromosome=chromosome,
        )
        m1 = Path(config.temporary_dir) / f"{STAGES['m1'].output_prefix}.{chromosome}.beta.gz"
        m2 = Path(config.temporary_dir) / f"{STAGES['m2'].output_prefix}.{chromosome}.beta.gz"
        m3 = Path(config.temporary_dir) / f"{stage.output_prefix}.{chromosome}.fit.tsv"
        score_prefix = Path(config.temporary_dir) / f"{stage.output_prefix}.grid.{chromosome}"
        score_path = Path(f"{score_prefix}.sscore")
        complete = checkpoint_matches(marker, config_hash) and m3.is_file() and score_path.is_file()
        if not (config.checkpoint.enabled and config.checkpoint.resume and complete):
            if not m1.is_file() or not m2.is_file():
                raise HerculesError(f"M3 requires M1/M2 posterior files: {m1}, {m2}")
            integrate_posterior_tables(
                m1,
                m2,
                m3,
                max_iter=int(config.m3.get("max_iter", 1000)),
                tol=float(config.m3.get("tol", 1e-6)),
            )
            score_path = run_plink2_score(
                plink2=config.tools.plink2,
                genotype_prefix=_format_path(
                    target_genotype,
                    chromosome=chromosome,
                    ancestry=config.target_ancestry,
                ),
                score_file=m3,
                output_prefix=score_prefix,
                grid=True,
            )
            _write_checkpoint(marker, config_hash, "m3", chromosome)
        grid_scores.append(score_path)

    aggregate_score_files(
        grid_scores,
        Path(config.output_dir) / f"{stage.output_prefix}.scores.tsv",
    )


def _run_ensemble(config: HerculesConfig) -> None:
    stage = STAGES["ensemble"]
    config_hash = configuration_hash(config.as_dict())
    marker = checkpoint_path(
        config.output_dir,
        stage=stage,
        trait=config.trait_name,
        chromosome=None,
    )
    output_prefix = Path(config.output_dir) / stage.output_prefix
    predictions = Path(f"{output_prefix}.predictions.tsv")
    metrics = Path(f"{output_prefix}.metrics.tsv")
    if (
        config.checkpoint.enabled
        and config.checkpoint.resume
        and checkpoint_matches(marker, config_hash)
        and predictions.is_file()
        and metrics.is_file()
    ):
        return

    sources = [
        ("M1_selected", Path(config.output_dir) / "HERCULES_M1.scores.tsv"),
        ("M2_selected", Path(config.output_dir) / "HERCULES_M2.scores.tsv"),
        ("M3", Path(config.output_dir) / "HERCULES_M3.scores.tsv"),
        ("M1_grid", Path(config.output_dir) / "HERCULES_M1.grid-scores.tsv"),
        ("M2_grid", Path(config.output_dir) / "HERCULES_M2.grid-scores.tsv"),
    ]
    combined: pd.DataFrame | None = None
    for label, path in sources:
        if not path.is_file():
            raise HerculesError(f"Ensemble input is missing: {path}")
        table = pd.read_csv(path, sep="\t")
        renamed = table.rename(
            columns={column: f"{label}_{column}" for column in table.columns if column != "IID"}
        )
        if combined is None:
            combined = renamed
        else:
            combined = combined.merge(renamed, on="IID", how="inner", sort=False)
    assert combined is not None
    predictor_path = Path(config.output_dir) / f"{stage.output_prefix}.inputs.tsv"
    combined.to_csv(predictor_path, sep="\t", index=False)

    result = run_process(
        [
            config.tools.rscript,
            ensemble_script_path(),
            predictor_path,
            config.inputs.phenotype_file,
            config.inputs.phenotype_column,
            ",".join(config.inputs.covariates),
            config.inputs.trait_type,
            output_prefix,
            str(config.execution.seed),
        ]
    )
    Path(f"{output_prefix}.stdout.log").write_text(result.stdout, encoding="utf-8")
    Path(f"{output_prefix}.stderr.log").write_text(result.stderr, encoding="utf-8")
    if not predictions.is_file() or not metrics.is_file():
        raise HerculesError("The R ensemble completed without producing predictions and metrics")
    _write_checkpoint(marker, config_hash, "ensemble", None)


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
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
