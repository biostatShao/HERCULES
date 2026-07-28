from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hercules.sumstats import prepare_fastgwa_sumstats, validate_fastgwa_header


def _table(*, include_prior: bool = True) -> pd.DataFrame:
    data: dict[str, object] = {
        "CHR": [22, 22],
        "SNP": ["rs1", "rs2"],
        "POS": [100, 200],
        "A1": ["A", "G"],
        "A2": ["C", "T"],
        "N": [1000, 1000],
        "AF1": [0.2, 0.3],
        "BETA": [0.01, -0.02],
        "SE": [0.03, 0.04],
        "P": [0.74, 0.62],
    }
    if include_prior:
        data["var_prior"] = [0.02, 0.04]
    return pd.DataFrame(data)


def test_var_prior_replaces_only_internal_p_column(tmp_path: Path) -> None:
    source = tmp_path / "input.fastGWA.tsv"
    destination = tmp_path / "internal.fastGWA.tsv"
    original = _table()
    original.to_csv(source, sep="\t", index=False)

    result = prepare_fastgwa_sumstats(source, destination)
    prepared = pd.read_csv(destination, sep="\t")

    assert result.used_var_prior is True
    assert result.variants == 2
    np.testing.assert_allclose(prepared["P"], [0.02, 0.04])
    np.testing.assert_allclose(pd.read_csv(source, sep="\t")["P"], [0.74, 0.62])
    assert "var_prior" not in prepared.columns


def test_missing_var_prior_uses_unit_variance(tmp_path: Path) -> None:
    source = tmp_path / "input.fastGWA.tsv"
    destination = tmp_path / "internal.fastGWA.tsv"
    _table(include_prior=False).to_csv(source, sep="\t", index=False)

    result = prepare_fastgwa_sumstats(source, destination)

    assert result.used_var_prior is False
    np.testing.assert_array_equal(pd.read_csv(destination, sep="\t")["P"], [1.0, 1.0])


@pytest.mark.parametrize("invalid", [0.0, -0.1, np.nan, np.inf])
def test_invalid_var_prior_is_rejected(tmp_path: Path, invalid: float) -> None:
    source = tmp_path / "invalid.fastGWA.tsv"
    table = _table()
    table.loc[1, "var_prior"] = invalid
    table.to_csv(source, sep="\t", index=False)

    with pytest.raises(ValueError, match="var_prior.*greater than zero"):
        prepare_fastgwa_sumstats(source, tmp_path / "internal.tsv")


def test_missing_required_column_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "missing.fastGWA.tsv"
    _table().drop(columns="SE").to_csv(source, sep="\t", index=False)

    with pytest.raises(ValueError, match="missing required columns: SE"):
        validate_fastgwa_header(source)


def test_association_p_value_keeps_normal_range_semantics(tmp_path: Path) -> None:
    source = tmp_path / "bad-p.fastGWA.tsv"
    table = _table()
    table.loc[0, "P"] = 1.2
    table.to_csv(source, sep="\t", index=False)

    with pytest.raises(ValueError, match="column 'P'.*between zero and one"):
        prepare_fastgwa_sumstats(source, tmp_path / "internal.tsv")
