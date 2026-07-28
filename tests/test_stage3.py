from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from hercules.exceptions import HerculesError
from hercules.resources import ensemble_script_path
from hercules.workflow import (
    _aggregate_single_score,
    _build_stage3_score_matrix,
    _ensure_disjoint_iids,
)


def test_stage3_predictor_matrix_has_exactly_two_genetic_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = SimpleNamespace(
        chromosomes=("22",),
        target_ancestry="TARGET",
        temporary_dir=str(tmp_path),
        tools=SimpleNamespace(plink2="plink2"),
    )
    (tmp_path / "HERCULES_M1.22.fit.gz").write_text("fixture", encoding="utf-8")
    (tmp_path / "HERCULES_M3.22.posterior.tsv").write_text("fixture", encoding="utf-8")
    observed_score_files: list[str] = []

    def fake_score(**kwargs):
        observed_score_files.append(Path(kwargs["score_file"]).name)
        output = Path(f'{kwargs["output_prefix"]}.sscore')
        output.write_text("fixture", encoding="utf-8")
        return output

    def fake_aggregate(paths, intermediate_path, output_column):
        return pd.DataFrame(
            {"IID": ["sample1", "sample2"], output_column: [0.1, 0.2]}
        )

    monkeypatch.setattr("hercules.workflow.run_plink2_score", fake_score)
    monkeypatch.setattr("hercules.workflow._aggregate_single_score", fake_aggregate)

    result = _build_stage3_score_matrix(
        config, cohort="validation", genotype="target_validation"
    )

    assert list(result.columns) == [
        "IID", "target_stage1_score", "calibrated_stage2_score"
    ]
    assert observed_score_files == [
        "HERCULES_M1.22.fit.gz",
        "HERCULES_M3.22.posterior.tsv",
    ]
    assert not any("M2" in name or "grid" in name.lower() for name in observed_score_files)


def test_validation_and_test_iids_must_be_disjoint() -> None:
    validation = pd.DataFrame(
        {
            "IID": ["v1", "shared"],
            "target_stage1_score": [0.1, 0.2],
            "calibrated_stage2_score": [0.3, 0.4],
        }
    )
    test = pd.DataFrame(
        {
            "IID": ["shared", "t1"],
            "target_stage1_score": [0.5, 0.6],
            "calibrated_stage2_score": [0.7, 0.8],
        }
    )

    with pytest.raises(HerculesError, match="overlapping IID detected: shared"):
        _ensure_disjoint_iids(validation, test)


def test_plink_dosage_diagnostic_is_not_a_stage3_predictor(tmp_path: Path) -> None:
    score = tmp_path / "score.sscore"
    score.write_text(
        "#FID\tIID\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_SUM\n"
        "sample1\tsample1\t12\t0.25\n",
        encoding="utf-8",
    )

    result = _aggregate_single_score(
        [score], tmp_path / "aggregate.tsv", "target_stage1_score"
    )

    assert list(result.columns) == ["IID", "target_stage1_score"]
    assert result.loc[0, "target_stage1_score"] == pytest.approx(0.25)


def test_chromosome_score_sample_sets_must_match(tmp_path: Path) -> None:
    first = tmp_path / "chr1.sscore"
    second = tmp_path / "chr2.sscore"
    first.write_text(
        "#FID\tIID\tSCORE1_SUM\n0\ta\t0.1\n0\tb\t0.2\n",
        encoding="utf-8",
    )
    second.write_text(
        "#FID\tIID\tSCORE1_SUM\n0\ta\t0.3\n0\tc\t0.4\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="sample sets differ across chromosomes"):
        _aggregate_single_score(
            [first, second], tmp_path / "aggregate.tsv", "target_stage1_score"
        )


def test_target_and_calibrated_score_sample_sets_must_match(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = SimpleNamespace(
        chromosomes=("22",),
        target_ancestry="TARGET",
        temporary_dir=str(tmp_path),
        tools=SimpleNamespace(plink2="plink2"),
    )
    (tmp_path / "HERCULES_M1.22.fit.gz").write_text("fixture", encoding="utf-8")
    (tmp_path / "HERCULES_M3.22.posterior.tsv").write_text("fixture", encoding="utf-8")
    monkeypatch.setattr(
        "hercules.workflow.run_plink2_score",
        lambda **kwargs: Path(f'{kwargs["output_prefix"]}.sscore'),
    )

    def mismatched_scores(paths, intermediate_path, output_column):
        iids = ["sample1", "sample2"] if output_column.startswith("target") else ["sample1"]
        return pd.DataFrame({"IID": iids, output_column: [0.1] * len(iids)})

    monkeypatch.setattr(
        "hercules.workflow._aggregate_single_score", mismatched_scores
    )

    with pytest.raises(HerculesError, match="score sample sets differ"):
        _build_stage3_score_matrix(
            config, cohort="validation", genotype="target_validation"
        )


def test_superlearner_libraries_and_test_outcome_boundary() -> None:
    script = ensemble_script_path().read_text(encoding="utf-8")

    assert 'c("SL.lasso.HERCULES", "SL.ridge.HERCULES", "SL.nnet")' in script
    assert 'c("SL.lasso.HERCULES", "SL.nnet")' in script
    assert 'meta_method <- "method.AUC"' in script
    assert "best_individual" not in script
    assert "SL.glmnet" in script and "alpha = 1" in script and "alpha = 0" in script
    assert "model$times <- NULL" in script
    assert "test_rows <- match(predictions$IID, test_phenotype$IID)" in script
    assert "validation_rows <- match(validation_scores$IID, validation_phenotype$IID)" in script
    assert "na.omit" not in script
    assert "covariate_model" not in script
    assert "residuals(" not in script
    assert "Test phenotype contains missing or non-finite values" in script
    assert "Binary test phenotype must contain both 0 and 1" in script
    assert script.index("model <- SuperLearner(") < script.index(
        "test_phenotype <- fread(test_phenotype_path)"
    )
