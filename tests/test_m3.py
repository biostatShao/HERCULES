from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hercules.m3 import integrate_posterior_tables, vi_bayes_paper


def test_vi_bayes_paper_is_deterministic_and_finite() -> None:
    mu_target = np.array([0.02, -0.01, 0.005, 0.03])
    variance_target = np.array([0.001, 0.002, 0.0015, 0.003])
    mu_base = np.array([0.015, -0.006, 0.008, 0.025])
    variance_base = np.array([0.0012, 0.0018, 0.0011, 0.0028])
    first = vi_bayes_paper(mu_target, variance_target, mu_base, variance_base, max_iter=1000)
    second = vi_bayes_paper(mu_target, variance_target, mu_base, variance_base, max_iter=1000)
    np.testing.assert_array_equal(first.beta_target, second.beta_target)
    assert np.isfinite(first.beta_target).all()
    assert np.isfinite(first.tau_squared)


def test_integrate_posterior_tables_preserves_100_candidate_contract(tmp_path: Path) -> None:
    keys = pd.DataFrame(
        {
            "CHR": [1, 1],
            "SNP": ["rs2", "rs1"],
            "POS": [20, 10],
            "A1": ["A", "C"],
            "A2": ["G", "T"],
        }
    )
    m1 = pd.concat(
        [
            keys,
            pd.DataFrame({f"BETA_{index}": [0.001 * (index + 1)] * 2 for index in range(100)}),
            pd.DataFrame({f"VAR_BETA_{index}": [0.01] * 2 for index in range(100)}),
        ],
        axis=1,
    )
    m2 = pd.concat(
        [
            keys,
            pd.DataFrame({f"BETA_{index}": [0.0015 * (index + 1)] * 2 for index in range(100)}),
            pd.DataFrame({f"VAR_BETA_{index}": [0.012] * 2 for index in range(100)}),
        ],
        axis=1,
    )
    m1_path = tmp_path / "m1.tsv"
    m2_path = tmp_path / "m2.tsv"
    output = tmp_path / "m3.tsv"
    m1.to_csv(m1_path, sep="\t", index=False)
    m2.to_csv(m2_path, sep="\t", index=False)
    integrate_posterior_tables(m1_path, m2_path, output)
    result = pd.read_csv(output, sep="\t")
    assert list(result.columns[:5]) == ["CHR", "SNP", "POS", "A1", "A2"]
    assert list(result["SNP"]) == ["rs1", "rs2"]
    assert [column for column in result if column.startswith("BETA_")] == [
        f"BETA_{index}" for index in range(100)
    ]


def test_vi_bayes_paper_rejects_zero_variance_inputs() -> None:
    values = np.array([0.01, 0.02])
    with pytest.raises(ValueError, match="strictly positive"):
        vi_bayes_paper(values, np.array([0.0, 0.01]), values, np.array([0.01, 0.01]))
