from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad

import hercules.m3 as m3
from hercules.m3 import (
    LAMBDA_PRIOR,
    M3Result,
    calibrate_directional,
    integrate_posterior_tables,
    mean_field_eta_parameters,
    stage2_marginal_log_likelihood,
)


def _scalar_reference(
    target_beta: float,
    target_variance: float,
    base_beta: float,
    base_variance: float,
    *,
    max_iter: int = 1000,
    tol: float = 1e-10,
) -> tuple[float, float, float]:
    mean = target_beta
    variance = target_variance
    previous = -np.inf
    lambda_mean = 0.5

    for _ in range(max_iter):
        second_moment = variance + mean**2

        def log_kernel(value: float) -> float:
            return (
                -np.log(value)
                - 0.5
                / base_variance
                * (
                    second_moment / value**2
                    - 2.0 * base_beta * mean / value
                    + base_beta**2
                )
            )

        grid = np.linspace(1e-5, 1.0, 10001)
        shift = float(
            np.max(
                -np.log(grid)
                - 0.5
                / base_variance
                * (
                    second_moment / grid**2
                    - 2.0 * base_beta * mean / grid
                    + base_beta**2
                )
            )
        )

        def scaled(value: float) -> float:
            return np.exp(log_kernel(value) - shift)

        normalizer = quad(scaled, 0.0, 1.0, epsabs=1e-12, limit=300)[0]
        inverse = (
            quad(lambda value: scaled(value) / value, 0.0, 1.0, epsabs=1e-12, limit=300)[0]
            / normalizer
        )
        inverse_squared = (
            quad(lambda value: scaled(value) / value**2, 0.0, 1.0, epsabs=1e-12, limit=300)[0]
            / normalizer
        )
        lambda_mean = (
            quad(lambda value: scaled(value) * value, 0.0, 1.0, epsabs=1e-12, limit=300)[0]
            / normalizer
        )
        precision = 1.0 / target_variance + inverse_squared / base_variance
        new_variance = 1.0 / precision
        new_mean = new_variance * (
            target_beta / target_variance + base_beta * inverse / base_variance
        )
        change = max(abs(new_mean - mean), abs(new_variance - variance))
        mean, variance = new_mean, new_variance
        if change <= tol and abs(change - previous) <= tol:
            break
        previous = change

    return mean, variance, lambda_mean


def _posterior(beta: list[float], variance: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CHR": [22] * len(beta),
            "SNP": [f"rs{index + 1}" for index in range(len(beta))],
            "POS": [100 * (index + 1) for index in range(len(beta))],
            "A1": ["A"] * len(beta),
            "A2": ["G"] * len(beta),
            "BETA": beta,
            "VAR_BETA": variance,
        }
    )


def test_single_snp_matches_independent_quadrature_reference() -> None:
    expected = _scalar_reference(0.08, 0.03, 0.12, 0.05)
    result = calibrate_directional(
        np.array([0.08]),
        np.array([0.03]),
        np.array([0.12]),
        np.array([0.05]),
        max_iter=1000,
        tol=1e-10,
        quadrature_points=96,
    )

    np.testing.assert_allclose(result.beta[0], expected[0], rtol=2e-5, atol=2e-7)
    np.testing.assert_allclose(result.variance[0], expected[1], rtol=2e-5, atol=2e-7)
    np.testing.assert_allclose(result.lambda_mean[0], expected[2], rtol=2e-5, atol=2e-7)


def test_integrated_stage2_likelihood_matches_methods_formula() -> None:
    target_beta = np.array([0.08, -0.03])
    target_variance = np.array([0.03, 0.04])
    base_beta = np.array([0.12, 0.05])
    base_variance = np.array([0.05, 0.06])
    lambda_value = np.array([0.25, 0.8])

    marginal_variance = target_variance + lambda_value**2 * base_variance
    expected = -0.5 * (
        np.log(2.0 * np.pi * marginal_variance)
        + (target_beta - lambda_value * base_beta) ** 2 / marginal_variance
    )

    np.testing.assert_allclose(
        stage2_marginal_log_likelihood(
            target_beta,
            target_variance,
            base_beta,
            base_variance,
            lambda_value,
        ),
        expected,
        rtol=0.0,
        atol=0.0,
    )


def test_eta_coordinate_update_matches_mean_field_formula() -> None:
    target_beta = np.array([0.08, -0.03])
    target_variance = np.array([0.03, 0.04])
    base_beta = np.array([0.12, 0.05])
    base_variance = np.array([0.05, 0.06])
    expected_inverse_lambda = np.array([2.0, 1.5])
    expected_inverse_lambda_squared = np.array([5.0, 2.5])

    expected_variance = 1.0 / (
        1.0 / target_variance
        + expected_inverse_lambda_squared / base_variance
    )
    expected_mean = expected_variance * (
        target_beta / target_variance
        + base_beta * expected_inverse_lambda / base_variance
    )

    mean, variance = mean_field_eta_parameters(
        target_beta,
        target_variance,
        base_beta,
        base_variance,
        expected_inverse_lambda,
        expected_inverse_lambda_squared,
    )
    np.testing.assert_allclose(mean, expected_mean, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(variance, expected_variance, rtol=0.0, atol=0.0)


@pytest.mark.parametrize(
    ("target_variance", "base_variance", "lambda_value", "message"),
    [
        (0.0, 0.05, 0.5, "posterior variances"),
        (0.03, np.inf, 0.5, "posterior variances"),
        (0.03, 0.05, -0.1, "lambda values"),
        (0.03, 0.05, 1.1, "lambda values"),
    ],
)
def test_integrated_stage2_likelihood_rejects_invalid_inputs(
    target_variance: float,
    base_variance: float,
    lambda_value: float,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        stage2_marginal_log_likelihood(
            np.array([0.08]),
            np.array([target_variance]),
            np.array([0.12]),
            np.array([base_variance]),
            np.array([lambda_value]),
        )


def test_multi_snp_result_is_finite_deterministic_and_converged() -> None:
    args = (
        np.array([0.05, -0.02, 0.11]),
        np.array([0.02, 0.03, 0.04]),
        np.array([0.08, 0.01, -0.04]),
        np.array([0.05, 0.06, 0.07]),
    )
    first = calibrate_directional(*args, max_iter=1000, tol=1e-8)
    second = calibrate_directional(*args, max_iter=1000, tol=1e-8)

    for name in ("beta", "variance", "lambda_mean", "elbo"):
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
        assert np.isfinite(getattr(first, name)).all()
    assert first.converged.all()
    assert np.all((first.lambda_mean >= 0.0) & (first.lambda_mean <= 1.0))
    assert np.all(first.iterations > 0)


def test_model_is_directional() -> None:
    target_beta = np.array([0.03, -0.07])
    target_variance = np.array([0.02, 0.04])
    base_beta = np.array([0.15, 0.02])
    base_variance = np.array([0.08, 0.03])

    forward = calibrate_directional(
        target_beta, target_variance, base_beta, base_variance
    )
    reverse = calibrate_directional(
        base_beta, base_variance, target_beta, target_variance
    )

    assert not np.allclose(forward.beta, reverse.beta)


def test_uniform_prior_is_fixed_scientific_constant() -> None:
    assert LAMBDA_PRIOR == "Uniform(0,1)"


def test_selected_tables_are_accepted_without_grid_columns(tmp_path: Path) -> None:
    target = _posterior([0.01, 0.02], [0.03, 0.04])
    base = _posterior([0.015, -0.01], [0.05, 0.06])
    target_path = tmp_path / "target.fit.gz"
    base_path = tmp_path / "base.fit.gz"
    output = tmp_path / "m3.tsv"
    diagnostics = tmp_path / "diagnostics.tsv"
    target.to_csv(target_path, sep="\t", index=False)
    base.to_csv(base_path, sep="\t", index=False)

    integrate_posterior_tables(
        target_path, base_path, output, diagnostics_path=diagnostics
    )

    posterior = pd.read_csv(output, sep="\t")
    diagnostic_table = pd.read_csv(diagnostics, sep="\t")
    assert list(posterior.columns) == [
        "CHR", "SNP", "POS", "A1", "A2", "BETA", "VAR_BETA"
    ]
    assert list(diagnostic_table.columns) == [
        "CHR", "SNP", "POS", "A1", "A2",
        "LAMBDA_MEAN", "CONVERGED", "ITERATIONS", "ELBO",
    ]


def test_var_beta_is_passed_as_variance_without_squaring(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = _posterior([0.01], [0.2])
    base = _posterior([0.02], [0.3])
    target_path = tmp_path / "target.tsv"
    base_path = tmp_path / "base.tsv"
    target.to_csv(target_path, sep="\t", index=False)
    base.to_csv(base_path, sep="\t", index=False)
    captured: dict[str, np.ndarray] = {}

    def fake_calibration(target_beta, target_variance, base_beta, base_variance, **kwargs):
        captured["target"] = np.asarray(target_variance)
        captured["base"] = np.asarray(base_variance)
        return M3Result(
            beta=np.asarray(target_beta),
            variance=np.asarray(target_variance),
            lambda_mean=np.array([0.5]),
            converged=np.array([True]),
            iterations=np.array([1]),
            elbo=np.array([0.0]),
        )

    monkeypatch.setattr(m3, "calibrate_directional", fake_calibration)
    integrate_posterior_tables(target_path, base_path, tmp_path / "output.tsv")

    np.testing.assert_array_equal(captured["target"], [0.2])
    np.testing.assert_array_equal(captured["base"], [0.3])


def test_output_changes_when_only_base_effect_changes() -> None:
    target_beta = np.array([0.05])
    target_variance = np.array([0.02])
    first = calibrate_directional(
        target_beta, target_variance, np.array([0.01]), np.array([0.04])
    )
    second = calibrate_directional(
        target_beta, target_variance, np.array([0.3]), np.array([0.04])
    )
    assert not np.isclose(first.beta[0], second.beta[0])


def test_pairwise_api_does_not_accept_an_additional_donor(tmp_path: Path) -> None:
    parameters = inspect.signature(integrate_posterior_tables).parameters
    positional = [
        parameter
        for parameter in parameters.values()
        if parameter.kind in {
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        }
    ]
    assert [parameter.name for parameter in positional] == [
        "target_path", "base_path", "output_path"
    ]


def test_unconverged_m3_is_not_written_as_a_successful_posterior(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target_path = tmp_path / "target.tsv"
    base_path = tmp_path / "base.tsv"
    output_path = tmp_path / "output.tsv"
    diagnostics_path = tmp_path / "diagnostics.tsv"
    _posterior([0.01], [0.03]).to_csv(target_path, sep="\t", index=False)
    _posterior([0.02], [0.04]).to_csv(base_path, sep="\t", index=False)

    monkeypatch.setattr(
        "hercules.m3.calibrate_directional",
        lambda *args, **kwargs: M3Result(
            beta=np.array([0.015]),
            variance=np.array([0.02]),
            lambda_mean=np.array([0.5]),
            converged=np.array([False]),
            iterations=np.array([1000]),
            elbo=np.array([-1.0]),
        ),
    )

    with pytest.raises(RuntimeError, match="did not converge"):
        integrate_posterior_tables(
            target_path,
            base_path,
            output_path,
            diagnostics_path=diagnostics_path,
        )

    assert diagnostics_path.is_file()
    assert not output_path.exists()


@pytest.mark.parametrize(
    "problem", ["duplicate", "identity", "alleles", "variance", "empty"]
)
def test_alignment_and_input_failures_are_actionable(tmp_path: Path, problem: str) -> None:
    target = _posterior([0.01, 0.02], [0.03, 0.04])
    base = _posterior([0.015, -0.01], [0.05, 0.06])
    expected = ""
    if problem == "duplicate":
        target.loc[1, ["SNP", "POS"]] = target.loc[0, ["SNP", "POS"]]
        expected = "duplicate SNP key"
    elif problem == "identity":
        base.loc[0, "POS"] += 1
        expected = "Incompatible target/base variant identity"
    elif problem == "alleles":
        base.loc[0, "A1"] = "T"
        expected = "Incompatible target/base alleles"
    elif problem == "variance":
        base.loc[0, "VAR_BETA"] = 0.0
        expected = "VAR_BETA must be finite and greater than zero"
    else:
        base["SNP"] = ["other1", "other2"]
        expected = "intersection is empty"

    target_path = tmp_path / "target.tsv"
    base_path = tmp_path / "base.tsv"
    target.to_csv(target_path, sep="\t", index=False)
    base.to_csv(base_path, sep="\t", index=False)

    with pytest.raises(ValueError, match=expected):
        integrate_posterior_tables(target_path, base_path, tmp_path / "output.tsv")
