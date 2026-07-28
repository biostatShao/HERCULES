from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from hercules.core.model.gridsearch.HerculesGrid import HerculesGrid
from hercules.core.model.per_snp_prior import prepare_initial_per_snp_precision
from hercules.workflow import _per_snp_prior_cli_args


def test_variance_is_converted_to_column_precision() -> None:
    precision = prepare_initial_per_snp_precision(
        np.array([0.02, 0.04, 0.08]),
        input_type="variance",
        column="PVAL",
        chromosome="22",
        float_precision="float32",
    )

    assert precision.shape == (3, 1)
    assert precision.dtype == np.float32
    np.testing.assert_allclose(
        precision[:, 0], np.array([50.0, 25.0, 12.5], dtype=np.float32)
    )


@pytest.mark.parametrize("invalid", [0.0, -0.1, np.nan, np.inf])
def test_invalid_prior_variance_is_rejected(invalid: float) -> None:
    with pytest.raises(ValueError, match="finite values greater than zero"):
        prepare_initial_per_snp_precision(
            np.array([0.02, invalid]),
            input_type="variance",
            column="PVAL",
            chromosome="22",
            float_precision="float64",
        )


def test_stage_initial_prior_defaults_to_fastgwa_pval_variance() -> None:
    assert _per_snp_prior_cli_args({}) == [
        "--initial-per-snp-prior-column",
        "PVAL",
        "--per-snp-prior-input-type",
        "variance",
    ]


def test_stage_initial_prior_configuration_is_forwarded() -> None:
    assert _per_snp_prior_cli_args(
        {
            "per_snp_prior": {
                "enabled": True,
                "column": "ANNOTATION_VARIANCE",
                "input_type": "variance",
            }
        }
    ) == [
        "--initial-per-snp-prior-column",
        "ANNOTATION_VARIANCE",
        "--per-snp-prior-input-type",
        "variance",
    ]


def test_m_step_updates_tau_beta_after_per_snp_initialization() -> None:
    model = object.__new__(HerculesGrid)
    model.fix_params = {}
    model.pi = np.array([0.1, 0.2], dtype=np.float64)
    model.gdl = SimpleNamespace(m=2)
    model.zeta = {
        "22": np.array([[0.02, 0.04], [0.03, 0.06]], dtype=np.float64)
    }
    model.float_precision = "float64"
    model.tau_beta = {"22": np.array([[50.0], [25.0]], dtype=np.float64)}

    model.update_tau_beta()

    np.testing.assert_allclose(model.tau_beta, [4.0, 4.0])
