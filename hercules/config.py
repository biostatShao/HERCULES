"""Typed YAML configuration with deterministic precedence and validation."""

from __future__ import annotations

import json
import os
import shutil
from collections.abc import Mapping, MutableMapping, Sequence
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .exceptions import ConfigurationError
from .paths import safe_path_component
from .sumstats import validate_fastgwa_header


@dataclass(slots=True)
class SummaryStatisticsConfig:
    base_path: str = ""
    target_path: str = ""
    base_columns: dict[str, str] = field(default_factory=dict)
    target_columns: dict[str, str] = field(default_factory=dict)


@dataclass(slots=True)
class InputConfig:
    summary_statistics: SummaryStatisticsConfig = field(default_factory=SummaryStatisticsConfig)
    functional_annotation: str = ""
    per_snp_heritability: str = ""
    ld_reference: dict[str, str] = field(default_factory=dict)
    genotype_prefixes: dict[str, str] = field(default_factory=dict)
    validation_genotype: str = ""
    phenotype_file: str = ""
    target_validation_genotype: str = ""
    target_validation_phenotype: str = ""
    target_test_genotype: str = ""
    target_test_phenotype: str = ""
    phenotype_column: str = ""
    covariates: tuple[str, ...] = ()
    trait_type: str = "quantitative"


@dataclass(slots=True)
class ToolConfig:
    plink: str = "plink"
    plink2: str = "plink2"
    rscript: str = "Rscript"


@dataclass(slots=True)
class ExecutionConfig:
    threads: int = 1
    parallel_jobs: int = 1
    seed: int = 7209


@dataclass(slots=True)
class CheckpointConfig:
    enabled: bool = True
    resume: bool = True


@dataclass(slots=True)
class LoggingConfig:
    level: str = "INFO"
    file: str = ""


@dataclass(slots=True)
class HerculesConfig:
    trait_name: str = ""
    chromosomes: tuple[str, ...] = tuple(str(i) for i in range(1, 23))
    base_ancestry: str = ""
    target_ancestry: str = ""
    inputs: InputConfig = field(default_factory=InputConfig)
    output_dir: str = "results"
    temporary_dir: str = "results/tmp"
    tools: ToolConfig = field(default_factory=ToolConfig)
    execution: ExecutionConfig = field(default_factory=ExecutionConfig)
    m1: dict[str, Any] = field(default_factory=dict)
    m2: dict[str, Any] = field(default_factory=dict)
    m3: dict[str, Any] = field(
        default_factory=lambda: {
            "model": "directional_pairwise_uniform_lambda",
            "lambda_prior": "uniform_0_1",
            "max_iter": 1000,
            "tol": 1e-6,
            "quadrature_points": 32,
        }
    )
    ensemble: dict[str, Any] = field(
        default_factory=lambda: {
            "quantitative_learners": ["lasso", "ridge", "neural_network"],
            "binary_learners": ["lasso", "neural_network"],
            "binary_method": "method.AUC",
        }
    )
    checkpoint: CheckpointConfig = field(default_factory=CheckpointConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def validate(self, *, check_paths: bool = True) -> tuple[str, ...]:
        """Validate inexpensive conditions and return non-fatal warnings."""

        errors: list[str] = []
        warnings: list[str] = []

        if not self.trait_name.strip():
            errors.append("trait_name is required")
        else:
            try:
                safe_path_component(self.trait_name, "trait_name")
            except ValueError as exc:
                errors.append(str(exc))
        if not self.base_ancestry.strip():
            errors.append("base_ancestry is required")
        if not self.target_ancestry.strip():
            errors.append("target_ancestry is required")
        if not self.chromosomes:
            errors.append("at least one chromosome is required")
        if self.execution.threads < 1:
            errors.append("execution.threads must be at least 1")
        if self.execution.parallel_jobs < 1:
            errors.append("execution.parallel_jobs must be at least 1")
        if self.inputs.trait_type not in {"quantitative", "binary"}:
            errors.append("inputs.trait_type must be 'quantitative' or 'binary'")
        if not self.inputs.phenotype_column.strip():
            errors.append("inputs.phenotype_column is required")

        required_files = (
            (
                "inputs.summary_statistics.base_path",
                self.inputs.summary_statistics.base_path,
                self.base_ancestry,
            ),
            (
                "inputs.summary_statistics.target_path",
                self.inputs.summary_statistics.target_path,
                self.target_ancestry,
            ),
            (
                "inputs.target_validation_phenotype",
                self.inputs.target_validation_phenotype,
                self.target_ancestry,
            ),
        )
        for label, value, ancestry in required_files:
            if not value:
                errors.append(f"{label} is required")
            elif check_paths and not _configured_file_exists(
                value, self.chromosomes, ancestry=ancestry
            ):
                errors.append(f"{label} does not exist: {value}")

        if check_paths:
            for label, value, ancestry in (
                (
                    "inputs.summary_statistics.base_path",
                    self.inputs.summary_statistics.base_path,
                    self.base_ancestry,
                ),
                (
                    "inputs.summary_statistics.target_path",
                    self.inputs.summary_statistics.target_path,
                    self.target_ancestry,
                ),
            ):
                for resolved in _configured_values(
                    value, self.chromosomes, ancestry=ancestry
                ):
                    if not Path(resolved).is_file():
                        continue
                    try:
                        validate_fastgwa_header(resolved)
                    except ValueError as exc:
                        errors.append(f"{label}: {exc}")

        for ancestry, path in (
            (self.base_ancestry, self.inputs.ld_reference.get("base", "")),
            (self.target_ancestry, self.inputs.ld_reference.get("target", "")),
        ):
            if not path:
                errors.append(f"inputs.ld_reference is missing the {ancestry or 'configured'} ancestry entry")
            elif check_paths and not _configured_path_exists(
                path, self.chromosomes, ancestry=ancestry
            ):
                errors.append(f"LD reference does not exist: {path}")

        if not self.inputs.validation_genotype:
            errors.append("inputs.validation_genotype is required")
        elif check_paths and not _configured_genotype_exists(
            self.inputs.validation_genotype,
            self.chromosomes,
            ancestry=self.target_ancestry,
        ):
            errors.append(
                "inputs.validation_genotype does not resolve to a PLINK .bed/.pgen prefix: "
                f"{self.inputs.validation_genotype}"
            )

        target_genotype = self.inputs.genotype_prefixes.get("target", "")
        if target_genotype:
            warnings.append(
                "inputs.genotype_prefixes.target is retained only as legacy metadata; "
                "Stage 3 uses explicit target validation and test genotypes."
            )

        for label, genotype in (
            ("inputs.target_validation_genotype", self.inputs.target_validation_genotype),
            ("inputs.target_test_genotype", self.inputs.target_test_genotype),
        ):
            if not genotype:
                errors.append(f"{label} is required")
            elif check_paths and not _configured_genotype_exists(
                genotype,
                self.chromosomes,
                ancestry=self.target_ancestry,
            ):
                errors.append(
                    f"{label} does not resolve to a PLINK .bed/.pgen prefix: {genotype}"
                )

        if self.inputs.target_test_phenotype and check_paths and not _configured_file_exists(
            self.inputs.target_test_phenotype,
            self.chromosomes,
            ancestry=self.target_ancestry,
        ):
            errors.append(
                "inputs.target_test_phenotype does not exist: "
                f"{self.inputs.target_test_phenotype}"
            )

        for stage_name, model_config, default_genotype, ancestry in (
            ("m1", self.m1, self.inputs.validation_genotype, self.target_ancestry),
            (
                "m2",
                self.m2,
                self.inputs.genotype_prefixes.get(
                    "base_validation", self.inputs.validation_genotype
                ),
                self.base_ancestry,
            ),
        ):
            if str(model_config.get("hyperparameter_search", "grid")).lower() not in {"grid", "gs"}:
                errors.append(f"{stage_name}.hyperparameter_search must be 'grid'")
            if int(model_config.get("pi_steps", 10)) != 10:
                errors.append(f"{stage_name}.pi_steps must be 10 for the published workflow")
            if int(model_config.get("sigma_epsilon_steps", 10)) != 10:
                errors.append(
                    f"{stage_name}.sigma_epsilon_steps must be 10 for the published workflow"
                )
            if str(model_config.get("sumstats_format", "fastgwa")).lower() != "fastgwa":
                errors.append(
                    f"{stage_name}.sumstats_format must be 'fastgwa' for the "
                    "validated HERCULES workflow"
                )
            prior = model_config.get("per_snp_prior", {})
            if not isinstance(prior, Mapping):
                errors.append(f"{stage_name}.per_snp_prior must be a mapping")
            else:
                allowed_prior_keys = {
                    "enabled",
                    "source",
                    "column",
                    "input_type",
                    "update_during_inference",
                }
                unknown_prior_keys = set(prior) - allowed_prior_keys
                if unknown_prior_keys:
                    errors.append(
                        f"{stage_name}.per_snp_prior has unknown keys: "
                        + ", ".join(sorted(str(key) for key in unknown_prior_keys))
                    )
                if not isinstance(prior.get("enabled", True), bool):
                    errors.append(f"{stage_name}.per_snp_prior.enabled must be boolean")
                elif not prior.get("enabled", True):
                    errors.append(
                        f"{stage_name}.per_snp_prior.enabled must be true for the "
                        "validated HERCULES M1/M2 workflow"
                    )
                if prior.get("source", "summary_statistics") != "summary_statistics":
                    errors.append(
                        f"{stage_name}.per_snp_prior.source must be 'summary_statistics'"
                    )
                column = prior.get("column", "PVAL")
                if not isinstance(column, str) or not column.strip():
                    errors.append(f"{stage_name}.per_snp_prior.column must be a non-empty string")
                if prior.get("input_type", "variance") not in {"variance", "precision"}:
                    errors.append(
                        f"{stage_name}.per_snp_prior.input_type must be 'variance' or 'precision'"
                    )
                if prior.get("update_during_inference", True) is not True:
                    errors.append(
                        f"{stage_name}.per_snp_prior.update_during_inference must be true"
                    )

            if str(model_config.get("grid_metric", "validation")) != "validation":
                continue
            validation_phenotype = str(model_config.get("validation_phenotype", ""))
            if not validation_phenotype:
                errors.append(f"{stage_name}.validation_phenotype is required")
            elif check_paths and not _configured_file_exists(
                validation_phenotype,
                self.chromosomes,
                ancestry=ancestry,
            ):
                errors.append(
                    f"{stage_name}.validation_phenotype does not exist: "
                    f"{validation_phenotype}"
                )
            validation_genotype = str(
                model_config.get("validation_genotype", default_genotype)
            )
            if not validation_genotype:
                errors.append(f"{stage_name}.validation_genotype is required")
            elif check_paths and not _configured_genotype_exists(
                validation_genotype,
                self.chromosomes,
                ancestry=ancestry,
            ):
                errors.append(
                    f"{stage_name}.validation_genotype does not resolve to a PLINK "
                    f".bed/.pgen prefix: {validation_genotype}"
                )
            validation_keep = str(model_config.get("validation_keep", ""))
            if validation_keep and check_paths and not _configured_file_exists(
                validation_keep,
                self.chromosomes,
                ancestry=ancestry,
            ):
                errors.append(
                    f"{stage_name}.validation_keep does not exist: {validation_keep}"
                )

        if self.m3.get("model", "") != "directional_pairwise_uniform_lambda":
            errors.append("m3.model must be 'directional_pairwise_uniform_lambda'")
        if self.m3.get("lambda_prior", "") != "uniform_0_1":
            errors.append("m3.lambda_prior must be 'uniform_0_1'")
        if int(self.m3.get("quadrature_points", 32)) < 8:
            errors.append("m3.quadrature_points must be at least 8")
        if list(self.ensemble.get("quantitative_learners", [])) != [
            "lasso",
            "ridge",
            "neural_network",
        ]:
            errors.append(
                "ensemble.quantitative_learners must be lasso, ridge, neural_network"
            )
        if list(self.ensemble.get("binary_learners", [])) != [
            "lasso",
            "neural_network",
        ]:
            errors.append("ensemble.binary_learners must be lasso, neural_network")
        if self.ensemble.get("binary_method") != "method.AUC":
            errors.append("ensemble.binary_method must be 'method.AUC'")

        if check_paths:
            for label, executable in (
                ("tools.plink", self.tools.plink),
                ("tools.plink2", self.tools.plink2),
                ("tools.rscript", self.tools.rscript),
            ):
                if not executable_available(executable):
                    errors.append(f"{label} is not executable or discoverable on PATH: {executable}")

            for label, destination in (
                ("output_dir", self.output_dir),
                ("temporary_dir", self.temporary_dir),
            ):
                if not _writable_destination(destination):
                    errors.append(f"{label} is not writable or creatable: {destination}")

        if self.inputs.functional_annotation:
            if check_paths and not _configured_file_exists(
                self.inputs.functional_annotation, self.chromosomes
            ):
                errors.append(
                    "inputs.functional_annotation does not exist: "
                    f"{self.inputs.functional_annotation}"
                )
            warnings.append(
                "Raw functional annotations are recorded for provenance only; provide their derived "
                "positive SNP-specific initialization variances in each FastGWA var_prior column."
            )
        else:
            warnings.append(
                "No raw functional annotation file is configured; Stage 1 still uses var_prior when "
                "that column is present in the ancestry-specific FastGWA input."
            )

        if (
            self.inputs.per_snp_heritability
            and check_paths
            and not _configured_file_exists(
                self.inputs.per_snp_heritability,
                self.chromosomes,
            )
        ):
            errors.append(
                "inputs.per_snp_heritability does not exist: "
                f"{self.inputs.per_snp_heritability}"
            )

        if self.execution.threads > 1:
            warnings.append(
                "Multithreaded numerical equivalence is not established because the inherited native "
                "E-step has a possible OpenMP write race."
            )

        if errors:
            raise ConfigurationError("Invalid HERCULES configuration:\n- " + "\n- ".join(errors))
        return tuple(warnings)


DEFAULT_CONFIG: dict[str, Any] = HerculesConfig().as_dict()

_OPEN_MAPPING_KEYS = {
    "base_columns",
    "target_columns",
    "ld_reference",
    "genotype_prefixes",
    "m1",
    "m2",
    "m3",
    "ensemble",
}


def load_config(
    path: str | os.PathLike[str],
    *,
    environ: Mapping[str, str] | None = None,
    cli_overrides: Mapping[str, Any] | None = None,
) -> HerculesConfig:
    """Load config using: CLI override > environment > YAML > package default."""

    config_path = Path(path)
    if not config_path.is_file():
        raise ConfigurationError(f"Configuration file does not exist: {config_path}")
    raw = _load_yaml(config_path)
    if not isinstance(raw, Mapping):
        raise ConfigurationError("The configuration root must be a mapping")

    merged = deepcopy(DEFAULT_CONFIG)
    _deep_merge(merged, raw)
    _apply_environment(merged, os.environ if environ is None else environ)
    for dotted_key, value in (cli_overrides or {}).items():
        _set_dotted(merged, dotted_key, value)
    _validate_schema(merged)
    return _from_mapping(merged)


def _load_yaml(path: Path) -> Any:
    try:
        import yaml
    except ImportError as exc:
        # JSON is a strict subset of YAML and provides a useful dependency-free path.
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            raise ConfigurationError(
                "PyYAML is required to read non-JSON YAML configuration files"
            ) from exc
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"Could not parse configuration {path}: {exc}") from exc


def _deep_merge(target: MutableMapping[str, Any], source: Mapping[str, Any]) -> None:
    if not target:
        target.update(deepcopy(dict(source)))
        return
    for key, value in source.items():
        if key not in target:
            raise ConfigurationError(f"Unknown configuration key: {key}")
        if key in _OPEN_MAPPING_KEYS and isinstance(target[key], MutableMapping) and isinstance(value, Mapping):
            target[key].update(deepcopy(dict(value)))
            continue
        if isinstance(target[key], MutableMapping) and isinstance(value, Mapping):
            _deep_merge(target[key], value)
        else:
            target[key] = value


def _apply_environment(config: MutableMapping[str, Any], environ: Mapping[str, str]) -> None:
    prefix = "HERCULES__"
    for name, raw_value in environ.items():
        if not name.startswith(prefix):
            continue
        dotted_key = name[len(prefix) :].lower().replace("__", ".")
        _set_dotted(config, dotted_key, _parse_scalar(raw_value))


def _set_dotted(config: MutableMapping[str, Any], dotted_key: str, value: Any) -> None:
    parts = dotted_key.split(".")
    current: MutableMapping[str, Any] = config
    for part in parts[:-1]:
        next_value = current.get(part)
        if not isinstance(next_value, MutableMapping):
            raise ConfigurationError(f"Unknown configuration override: {dotted_key}")
        current = next_value
    if parts[-1] not in current and (len(parts) < 2 or parts[-2] not in _OPEN_MAPPING_KEYS):
        raise ConfigurationError(f"Unknown configuration override: {dotted_key}")
    current[parts[-1]] = value


def _parse_scalar(value: str) -> Any:
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def _from_mapping(data: Mapping[str, Any]) -> HerculesConfig:
    inputs = data["inputs"]
    summary = inputs["summary_statistics"]
    return HerculesConfig(
        trait_name=str(data["trait_name"]),
        chromosomes=tuple(str(c) for c in _as_sequence(data["chromosomes"])),
        base_ancestry=str(data["base_ancestry"]),
        target_ancestry=str(data["target_ancestry"]),
        inputs=InputConfig(
            summary_statistics=SummaryStatisticsConfig(
                base_path=str(summary["base_path"]),
                target_path=str(summary["target_path"]),
                base_columns=dict(summary["base_columns"]),
                target_columns=dict(summary["target_columns"]),
            ),
            functional_annotation=str(inputs["functional_annotation"]),
            per_snp_heritability=str(inputs["per_snp_heritability"]),
            ld_reference=dict(inputs["ld_reference"]),
            genotype_prefixes=dict(inputs["genotype_prefixes"]),
            validation_genotype=str(inputs["validation_genotype"]),
            phenotype_file=str(inputs["phenotype_file"]),
            target_validation_genotype=str(inputs["target_validation_genotype"]),
            target_validation_phenotype=str(inputs["target_validation_phenotype"]),
            target_test_genotype=str(inputs["target_test_genotype"]),
            target_test_phenotype=str(inputs["target_test_phenotype"]),
            phenotype_column=str(inputs["phenotype_column"]),
            covariates=tuple(str(c) for c in _as_sequence(inputs["covariates"])),
            trait_type=str(inputs["trait_type"]).lower(),
        ),
        output_dir=str(data["output_dir"]),
        temporary_dir=str(data["temporary_dir"]),
        tools=ToolConfig(**data["tools"]),
        execution=ExecutionConfig(**data["execution"]),
        m1=dict(data["m1"]),
        m2=dict(data["m2"]),
        m3=dict(data["m3"]),
        ensemble=dict(data["ensemble"]),
        checkpoint=CheckpointConfig(**data["checkpoint"]),
        logging=LoggingConfig(**data["logging"]),
    )


def _validate_schema(data: Mapping[str, Any]) -> None:
    for key in (
        "trait_name",
        "base_ancestry",
        "target_ancestry",
        "output_dir",
        "temporary_dir",
    ):
        _require_type(data[key], str, key)

    chromosomes = _require_sequence(data["chromosomes"], "chromosomes")
    for index, chromosome in enumerate(chromosomes):
        if isinstance(chromosome, bool) or not isinstance(chromosome, (str, int)):
            raise ConfigurationError(
                f"chromosomes[{index}] must be a string or integer, not {type(chromosome).__name__}"
            )

    inputs = _require_mapping(data["inputs"], "inputs")
    summary = _require_mapping(inputs["summary_statistics"], "inputs.summary_statistics")
    for key in ("base_path", "target_path"):
        _require_type(summary[key], str, f"inputs.summary_statistics.{key}")
    _require_string_mapping(summary["base_columns"], "inputs.summary_statistics.base_columns")
    _require_string_mapping(summary["target_columns"], "inputs.summary_statistics.target_columns")
    for key in (
        "functional_annotation",
        "per_snp_heritability",
        "validation_genotype",
        "phenotype_file",
        "target_validation_genotype",
        "target_validation_phenotype",
        "target_test_genotype",
        "target_test_phenotype",
        "phenotype_column",
        "trait_type",
    ):
        _require_type(inputs[key], str, f"inputs.{key}")
    _require_string_mapping(inputs["ld_reference"], "inputs.ld_reference")
    _require_string_mapping(inputs["genotype_prefixes"], "inputs.genotype_prefixes")
    for index, covariate in enumerate(_require_sequence(inputs["covariates"], "inputs.covariates")):
        _require_type(covariate, str, f"inputs.covariates[{index}]")

    tools = _require_mapping(data["tools"], "tools")
    for key in ("plink", "plink2", "rscript"):
        _require_type(tools[key], str, f"tools.{key}")

    execution = _require_mapping(data["execution"], "execution")
    for key in ("threads", "parallel_jobs", "seed"):
        value = execution[key]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ConfigurationError(f"execution.{key} must be an integer")

    for key in ("m1", "m2", "m3", "ensemble"):
        model_section = _require_mapping(data[key], key)
        _require_string_keys(model_section, key)
        if key in {"m1", "m2"} and "per_snp_prior" in model_section:
            prior_section = _require_mapping(
                model_section["per_snp_prior"], f"{key}.per_snp_prior"
            )
            _require_string_keys(prior_section, f"{key}.per_snp_prior")

    for section_name, keys in (("checkpoint", ("enabled", "resume")),):
        section = _require_mapping(data[section_name], section_name)
        for key in keys:
            _require_type(section[key], bool, f"{section_name}.{key}")

    logging = _require_mapping(data["logging"], "logging")
    for key in ("level", "file"):
        _require_type(logging[key], str, f"logging.{key}")


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{label} must be a mapping, not {type(value).__name__}")
    return value


def _require_sequence(value: Any, label: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{label} must be a sequence, not {type(value).__name__}")
    return value


def _require_type(value: Any, expected: type, label: str) -> None:
    if not isinstance(value, expected):
        raise ConfigurationError(f"{label} must be {expected.__name__}, not {type(value).__name__}")


def _require_string_mapping(value: Any, label: str) -> None:
    mapping = _require_mapping(value, label)
    for key, item in mapping.items():
        if not isinstance(key, str) or not isinstance(item, str):
            raise ConfigurationError(f"{label} keys and values must all be strings")


def _require_string_keys(value: Any, label: str) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ConfigurationError(f"{label} mapping keys must all be strings")
            _require_string_keys(item, f"{label}.{key}")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for index, item in enumerate(value):
            _require_string_keys(item, f"{label}[{index}]")


def _as_sequence(value: Any) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"Expected a sequence, received {type(value).__name__}")
    return value


def _genotype_prefix_exists(value: str) -> bool:
    path = Path(value)
    suffix = path.suffix.lower()
    if suffix in {".bed", ".pgen"}:
        prefix = path.with_suffix("")
    elif path.is_file():
        return False
    else:
        prefix = path
    bed_complete = all(Path(f"{prefix}{item}").is_file() for item in (".bed", ".bim", ".fam"))
    pgen_complete = all(Path(f"{prefix}{item}").is_file() for item in (".pgen", ".pvar", ".psam"))
    return bed_complete or pgen_complete


def _configured_values(
    value: str,
    chromosomes: Sequence[str],
    *,
    ancestry: str = "",
) -> tuple[str, ...]:
    if "{" not in value:
        return (value,)
    expanded: list[str] = []
    for chromosome in chromosomes:
        try:
            expanded.append(
                value.format(
                    chrom=chromosome,
                    chromosome=chromosome,
                    ancestry=ancestry,
                )
            )
        except KeyError:
            return (value,)
    return tuple(expanded)


def _configured_file_exists(
    value: str,
    chromosomes: Sequence[str],
    *,
    ancestry: str = "",
) -> bool:
    return all(
        Path(item).is_file()
        for item in _configured_values(value, chromosomes, ancestry=ancestry)
    )


def _configured_path_exists(
    value: str,
    chromosomes: Sequence[str],
    *,
    ancestry: str = "",
) -> bool:
    return all(
        Path(item).exists()
        for item in _configured_values(value, chromosomes, ancestry=ancestry)
    )


def _configured_genotype_exists(
    value: str,
    chromosomes: Sequence[str],
    *,
    ancestry: str = "",
) -> bool:
    return all(
        _genotype_prefix_exists(item)
        for item in _configured_values(value, chromosomes, ancestry=ancestry)
    )


def executable_available(value: str) -> bool:
    candidate = Path(value)
    if candidate.parent != Path(".") or candidate.is_absolute():
        return candidate.is_file() and os.access(candidate, os.X_OK)
    return shutil.which(value) is not None


def _writable_destination(value: str) -> bool:
    candidate = Path(value).expanduser()
    if candidate.exists():
        return candidate.is_dir() and os.access(candidate, os.W_OK)
    probe = candidate.parent
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return probe.is_dir() and os.access(probe, os.W_OK)
