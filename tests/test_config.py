from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from hercules.config import load_config
from hercules.exceptions import ConfigurationError


def test_load_and_validate_paths_with_spaces(config_path: Path) -> None:
    config = load_config(config_path)
    warnings = config.validate(check_paths=True)
    assert config.trait_name == "height"
    assert config.execution.seed == 7209
    assert any("annotation" in warning.lower() for warning in warnings)


def test_precedence_cli_over_environment_over_yaml(config_path: Path) -> None:
    config = load_config(
        config_path,
        environ={"HERCULES__EXECUTION__THREADS": "3", "HERCULES__EXECUTION__SEED": "99"},
        cli_overrides={"execution.threads": 5},
    )
    assert config.execution.threads == 5
    assert config.execution.seed == 99


def test_unknown_configuration_key_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text(json.dumps({"unexpected": True}), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="Unknown configuration key"):
        load_config(path)


def test_missing_required_value_is_actionable(config_path: Path) -> None:
    config = load_config(config_path, cli_overrides={"execution.threads": 0})
    with pytest.raises(ConfigurationError, match="execution.threads"):
        config.validate(check_paths=False)


def test_model_sections_accept_method_specific_parameters(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    config_mapping["m3"] = {"max_iter": 1000, "tol": 1e-6, "additional_rho": 0.25}
    path = tmp_path / "model.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    config = load_config(path)
    assert config.m3["additional_rho"] == 0.25


def test_m1_and_m2_accept_explicit_fixed_per_snp_prior(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    prior = {
        "enabled": True,
        "source": "summary_statistics",
        "column": "PVAL",
        "input_type": "variance",
        "fixed_during_inference": True,
    }
    config_mapping["m1"] = {**config_mapping["m1"], "per_snp_prior": prior}
    config_mapping["m2"] = {**config_mapping["m2"], "per_snp_prior": prior}
    path = tmp_path / "fixed-prior.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")

    config = load_config(path)

    assert config.m1["per_snp_prior"]["column"] == "PVAL"
    assert config.m2["per_snp_prior"]["fixed_during_inference"] is True


def test_non_fixed_per_snp_prior_is_rejected(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    config_mapping["m1"] = {
        **config_mapping["m1"],
        "per_snp_prior": {"fixed_during_inference": False},
    }
    path = tmp_path / "non-fixed-prior.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")

    with pytest.raises(ConfigurationError, match="fixed_during_inference must be true"):
        load_config(path).validate(check_paths=False)


@pytest.mark.parametrize(
    ("key", "value", "message"),
    [
        ("trait_name", None, "trait_name must be str"),
        ("execution", {"threads": "many"}, "execution.threads must be an integer"),
        ("chromosomes", "1,2", "chromosomes must be a sequence"),
    ],
)
def test_malformed_types_raise_configuration_errors(
    tmp_path: Path,
    config_mapping: dict[str, object],
    key: str,
    value: object,
    message: str,
) -> None:
    if key == "execution":
        config_mapping[key] = value
    else:
        config_mapping[key] = value
    path = tmp_path / "malformed.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match=message):
        load_config(path)


def test_environment_can_extend_model_specific_mapping(config_path: Path) -> None:
    config = load_config(config_path, environ={"HERCULES__M3__ADDITIONAL_RHO": "0.5"})
    assert config.m3["additional_rho"] == 0.5


def test_invalid_tool_path_is_rejected(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    tools = dict(config_mapping["tools"])
    tools["plink2"] = str(tmp_path / "missing-plink2")
    config_mapping["tools"] = tools
    path = tmp_path / "bad-tool.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="tools.plink2"):
        load_config(path).validate(check_paths=True)


def test_invalid_plink1_path_is_rejected(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    tools = dict(config_mapping["tools"])
    tools["plink"] = str(tmp_path / "missing-plink")
    config_mapping["tools"] = tools
    path = tmp_path / "bad-plink1.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="tools.plink"):
        load_config(path).validate(check_paths=True)


def test_configured_per_snp_file_is_validated(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    inputs = dict(config_mapping["inputs"])
    inputs["per_snp_heritability"] = str(tmp_path / "missing-per-snp.tsv")
    config_mapping["inputs"] = inputs
    path = tmp_path / "missing-per-snp.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="inputs.per_snp_heritability"):
        load_config(path).validate(check_paths=True)


def test_ancestry_placeholder_uses_stage_context(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    template = tmp_path / "sumstats_{ancestry}_chr_{chrom}.tsv"
    for ancestry in ("EUR", "AFR"):
        for chromosome in ("1", "2"):
            Path(str(template).format(ancestry=ancestry, chrom=chromosome)).write_text(
                "fixture\n", encoding="utf-8"
            )
    inputs = dict(config_mapping["inputs"])
    summary = dict(inputs["summary_statistics"])
    summary["base_path"] = str(template)
    summary["target_path"] = str(template)
    inputs["summary_statistics"] = summary
    config_mapping["inputs"] = inputs
    path = tmp_path / "ancestry-placeholder.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    load_config(path).validate(check_paths=True)


def test_trait_path_traversal_is_rejected(config_path: Path) -> None:
    config = load_config(config_path)
    config.trait_name = "../height"
    with pytest.raises(ConfigurationError, match="Unsafe trait_name"):
        config.validate(check_paths=False)


def test_incomplete_plink_fileset_is_rejected(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    incomplete = tmp_path / "incomplete"
    Path(f"{incomplete}.bed").write_bytes(b"")
    inputs = dict(config_mapping["inputs"])
    inputs["validation_genotype"] = str(incomplete)
    config_mapping["inputs"] = inputs
    path = tmp_path / "incomplete.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="PLINK .bed/.pgen prefix"):
        load_config(path).validate(check_paths=True)


def test_output_path_cannot_be_existing_file(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    output_file = tmp_path / "not-a-directory"
    output_file.write_text("file", encoding="utf-8")
    config_mapping["output_dir"] = str(output_file)
    path = tmp_path / "output-file.yaml"
    path.write_text(json.dumps(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="output_dir is not writable or creatable"):
        load_config(path).validate(check_paths=True)


def test_model_mapping_keys_must_be_strings(
    tmp_path: Path, config_mapping: dict[str, object]
) -> None:
    config_mapping["m1"] = {1: "invalid"}
    path = tmp_path / "non-string-key.yaml"
    path.write_text(yaml.safe_dump(config_mapping), encoding="utf-8")
    with pytest.raises(ConfigurationError, match="m1 mapping keys"):
        load_config(path)
