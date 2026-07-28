"""Directional pairwise Stage-2 calibration for HERCULES.

For each aligned SNP, the target effect is calibrated from one designated
base population under

``eta | lambda ~ Normal(lambda * b_base, lambda**2 * V_base)``,
``b_target | eta ~ Normal(eta, V_target)``, and
``lambda ~ Uniform(0, 1)``.

The mean-field family is ``q(eta) q(lambda)``. ``q(eta)`` is Gaussian and
``q(lambda)`` is a bounded density represented by deterministic Gauss-Legendre
quadrature on (0, 1). Coordinate updates maximize the directly evaluated ELBO.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp


M3_SCIENTIFIC_MODEL = "directional-pairwise-uniform-lambda-v1"
LAMBDA_PRIOR = "Uniform(0,1)"
VARIANT_KEYS = ("CHR", "SNP", "POS", "A1", "A2")
IDENTITY_KEYS = ("CHR", "SNP", "POS")
POSTERIOR_COLUMNS = (*VARIANT_KEYS, "BETA", "VAR_BETA")


@dataclass(frozen=True, slots=True)
class M3Result:
    beta: np.ndarray
    variance: np.ndarray
    lambda_mean: np.ndarray
    converged: np.ndarray
    iterations: np.ndarray
    elbo: np.ndarray


def stage2_marginal_log_likelihood(
    target_beta: np.ndarray,
    target_variance: np.ndarray,
    base_beta: np.ndarray,
    base_variance: np.ndarray,
    lambda_value: np.ndarray,
) -> np.ndarray:
    """Evaluate the Methods marginal Stage-2 likelihood for fixed lambda.

    Integrating out ``eta`` gives
    ``b_target | lambda ~ Normal(lambda*b_base,
    V_target + lambda**2*V_base)``.
    """

    try:
        (
            target_beta,
            target_variance,
            base_beta,
            base_variance,
            lambda_value,
        ) = np.broadcast_arrays(
            np.asarray(target_beta, dtype=np.float64),
            np.asarray(target_variance, dtype=np.float64),
            np.asarray(base_beta, dtype=np.float64),
            np.asarray(base_variance, dtype=np.float64),
            np.asarray(lambda_value, dtype=np.float64),
        )
    except ValueError as exc:
        raise ValueError("Stage-2 marginal-likelihood inputs are not broadcast-compatible") from exc
    if not np.isfinite(target_beta).all() or not np.isfinite(base_beta).all():
        raise ValueError("Stage-2 posterior means must be finite")
    if (
        not np.isfinite(target_variance).all()
        or not np.isfinite(base_variance).all()
        or np.any(target_variance <= 0.0)
        or np.any(base_variance <= 0.0)
    ):
        raise ValueError("Stage-2 posterior variances must be finite and greater than zero")
    if (
        not np.isfinite(lambda_value).all()
        or np.any(lambda_value < 0.0)
        or np.any(lambda_value > 1.0)
    ):
        raise ValueError("Stage-2 lambda values must be finite and within [0, 1]")
    marginal_variance = target_variance + np.square(lambda_value) * base_variance
    return -0.5 * (
        np.log(2.0 * np.pi * marginal_variance)
        + np.square(target_beta - lambda_value * base_beta) / marginal_variance
    )


def mean_field_eta_parameters(
    target_beta: np.ndarray,
    target_variance: np.ndarray,
    base_beta: np.ndarray,
    base_variance: np.ndarray,
    expected_inverse_lambda: np.ndarray,
    expected_inverse_lambda_squared: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the Gaussian ``q(eta)`` coordinate update from the Methods model."""

    precision = (
        1.0 / target_variance
        + expected_inverse_lambda_squared / base_variance
    )
    variance = 1.0 / precision
    mean = variance * (
        target_beta / target_variance
        + base_beta * expected_inverse_lambda / base_variance
    )
    return mean, variance


def calibrate_directional(
    target_beta: np.ndarray,
    target_variance: np.ndarray,
    base_beta: np.ndarray,
    base_variance: np.ndarray,
    *,
    max_iter: int = 1000,
    tol: float = 1e-6,
    quadrature_points: int = 32,
) -> M3Result:
    """Fit the directional target-from-base mean-field model.

    The Uniform(0, 1) prior has constant log density zero. The lambda density
    is optimized non-parametrically at fixed Gauss-Legendre nodes. Nodes do not
    include the singular endpoints, log-sum-exp normalization prevents
    underflow, and all input variances must be finite and strictly positive.
    """

    target_beta = _one_dimensional(target_beta, "target BETA")
    target_variance = _one_dimensional(target_variance, "target VAR_BETA")
    base_beta = _one_dimensional(base_beta, "base BETA")
    base_variance = _one_dimensional(base_variance, "base VAR_BETA")
    if not (
        target_beta.shape
        == target_variance.shape
        == base_beta.shape
        == base_variance.shape
    ):
        raise ValueError("Stage-2 posterior vectors must have identical shapes")
    if target_beta.size == 0:
        raise ValueError("Stage-2 posterior vectors must not be empty")
    if not np.isfinite(target_beta).all() or not np.isfinite(base_beta).all():
        raise ValueError("Stage-2 posterior means must be finite")
    for label, values in (
        ("target VAR_BETA", target_variance),
        ("base VAR_BETA", base_variance),
    ):
        if not np.isfinite(values).all() or np.any(values <= 0.0):
            raise ValueError(f"{label} must contain finite values greater than zero")
    if max_iter < 1:
        raise ValueError("max_iter must be at least one")
    if tol <= 0.0 or not np.isfinite(tol):
        raise ValueError("tol must be finite and greater than zero")
    if quadrature_points < 8:
        raise ValueError("quadrature_points must be at least eight")

    raw_nodes, raw_weights = leggauss(quadrature_points)
    lambda_nodes = 0.5 * (raw_nodes + 1.0)
    lambda_weights = 0.5 * raw_weights
    log_weights = np.log(lambda_weights)
    inverse_lambda = 1.0 / lambda_nodes
    inverse_lambda_squared = np.square(inverse_lambda)

    mean = target_beta.copy()
    variance = target_variance.copy()
    lambda_mean = np.full(target_beta.shape, 0.5, dtype=np.float64)
    elbo = np.full(target_beta.shape, -np.inf, dtype=np.float64)
    converged = np.zeros(target_beta.shape, dtype=bool)
    iterations = np.zeros(target_beta.shape, dtype=np.int32)
    active = np.ones(target_beta.shape, dtype=bool)

    log_two_pi = np.log(2.0 * np.pi)
    for iteration in range(1, max_iter + 1):
        second_moment = variance + np.square(mean)
        log_lambda_density_kernel = (
            -np.log(lambda_nodes)[None, :]
            - 0.5
            / base_variance[:, None]
            * (
                second_moment[:, None] * inverse_lambda_squared[None, :]
                - 2.0
                * base_beta[:, None]
                * mean[:, None]
                * inverse_lambda[None, :]
                + np.square(base_beta)[:, None]
            )
        )
        log_normalizer = logsumexp(
            log_lambda_density_kernel + log_weights[None, :], axis=1
        )
        log_mass = (
            log_lambda_density_kernel
            + log_weights[None, :]
            - log_normalizer[:, None]
        )
        mass = np.exp(log_mass)
        expected_inverse = mass @ inverse_lambda
        expected_inverse_squared = mass @ inverse_lambda_squared
        candidate_lambda_mean = mass @ lambda_nodes

        candidate_mean, candidate_variance = mean_field_eta_parameters(
            target_beta,
            target_variance,
            base_beta,
            base_variance,
            expected_inverse,
            expected_inverse_squared,
        )

        candidate_second_moment = candidate_variance + np.square(candidate_mean)
        expected_log_lambda = mass @ np.log(lambda_nodes)
        target_term = (
            -0.5 * (log_two_pi + np.log(target_variance))
            - 0.5
            * (np.square(target_beta - candidate_mean) + candidate_variance)
            / target_variance
        )
        base_term = (
            -0.5 * (log_two_pi + np.log(base_variance))
            - expected_log_lambda
            - 0.5
            / base_variance
            * (
                candidate_second_moment * expected_inverse_squared
                - 2.0 * base_beta * candidate_mean * expected_inverse
                + np.square(base_beta)
            )
        )
        eta_entropy = 0.5 * np.log(2.0 * np.pi * np.e * candidate_variance)
        log_density = log_lambda_density_kernel - log_normalizer[:, None]
        lambda_entropy = -np.sum(mass * log_density, axis=1)
        candidate_elbo = target_term + base_term + eta_entropy + lambda_entropy

        previously_active = active.copy()
        mean[previously_active] = candidate_mean[previously_active]
        variance[previously_active] = candidate_variance[previously_active]
        lambda_mean[previously_active] = candidate_lambda_mean[previously_active]
        iterations[previously_active] = iteration
        if iteration > 1:
            newly_converged = previously_active & (
                np.abs(candidate_elbo - elbo) <= tol
            )
            converged[newly_converged] = True
            active[newly_converged] = False
        elbo[previously_active] = candidate_elbo[previously_active]
        if not np.any(active):
            break

    return M3Result(
        beta=mean,
        variance=variance,
        lambda_mean=lambda_mean,
        converged=converged,
        iterations=iterations,
        elbo=elbo,
    )


def integrate_posterior_tables(
    target_path: str | Path,
    base_path: str | Path,
    output_path: str | Path,
    *,
    diagnostics_path: str | Path | None = None,
    max_iter: int = 1000,
    tol: float = 1e-6,
    quadrature_points: int = 32,
) -> Path:
    """Calibrate one selected target posterior from one selected base posterior."""

    target = _read_selected_posterior(target_path, "target")
    base = _read_selected_posterior(base_path, "base")
    shared_snps = target.merge(
        base,
        on="SNP",
        how="inner",
        suffixes=("_target", "_base"),
        sort=False,
        validate="one_to_one",
    )
    if not shared_snps.empty:
        identity_conflict = (
            (shared_snps["CHR_target"] != shared_snps["CHR_base"])
            | (shared_snps["POS_target"] != shared_snps["POS_base"])
        )
        if identity_conflict.any():
            row = shared_snps.loc[identity_conflict].iloc[0]
            raise ValueError(
                "Incompatible target/base variant identity for "
                f"{row['SNP']}: target={row['CHR_target']}:{row['POS_target']}, "
                f"base={row['CHR_base']}:{row['POS_base']}"
            )
    target = target.assign(_target_order=np.arange(len(target), dtype=np.int64))
    merged = target.merge(
        base,
        on=list(IDENTITY_KEYS),
        how="inner",
        suffixes=("_target", "_base"),
        sort=False,
        validate="one_to_one",
    )
    if merged.empty:
        raise ValueError(
            f"Selected target/base posterior intersection is empty: {target_path}, {base_path}"
        )
    incompatible = (
        (merged["A1_target"] != merged["A1_base"])
        | (merged["A2_target"] != merged["A2_base"])
    )
    if incompatible.any():
        row = merged.loc[incompatible].iloc[0]
        raise ValueError(
            "Incompatible target/base alleles for "
            f"{row['SNP']} at {row['CHR']}:{row['POS']}: "
            f"target={row['A1_target']}/{row['A2_target']}, "
            f"base={row['A1_base']}/{row['A2_base']}"
        )
    merged = merged.sort_values("_target_order", kind="stable").reset_index(drop=True)

    result = calibrate_directional(
        merged["BETA_target"].to_numpy(),
        merged["VAR_BETA_target"].to_numpy(),
        merged["BETA_base"].to_numpy(),
        merged["VAR_BETA_base"].to_numpy(),
        max_iter=max_iter,
        tol=tol,
        quadrature_points=quadrature_points,
    )
    keys = pd.DataFrame(
        {
            "CHR": merged["CHR"],
            "SNP": merged["SNP"],
            "POS": merged["POS"],
            "A1": merged["A1_target"],
            "A2": merged["A2_target"],
        }
    )
    posterior = keys.assign(BETA=result.beta, VAR_BETA=result.variance)
    diagnostics = keys.assign(
        LAMBDA_MEAN=result.lambda_mean,
        CONVERGED=result.converged,
        ITERATIONS=result.iterations,
        ELBO=result.elbo,
    )

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if diagnostics_path is not None:
        diagnostic_destination = Path(diagnostics_path)
        diagnostic_destination.parent.mkdir(parents=True, exist_ok=True)
        diagnostics.to_csv(diagnostic_destination, sep="\t", index=False)
    if not result.converged.all():
        failed = np.flatnonzero(~result.converged)
        example = merged.iloc[int(failed[0])]
        raise RuntimeError(
            "Directional M3 did not converge for "
            f"{len(failed)} variant(s); first failure is {example['SNP']} at "
            f"{example['CHR']}:{example['POS']}. Increase m3.max_iter or inspect "
            "the convergence diagnostics before retrying."
        )
    posterior.to_csv(destination, sep="\t", index=False)
    return destination


def _read_selected_posterior(path: str | Path, label: str) -> pd.DataFrame:
    source = Path(path)
    try:
        table = pd.read_csv(source, sep="\t", usecols=list(POSTERIOR_COLUMNS))
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise ValueError(f"Could not read selected {label} posterior {source}: {exc}") from exc
    if table.empty:
        raise ValueError(f"Selected {label} posterior is empty: {source}")
    if table["SNP"].duplicated().any() or table.duplicated(list(IDENTITY_KEYS)).any():
        duplicate = table.loc[
            table["SNP"].duplicated(keep=False) | table.duplicated(list(IDENTITY_KEYS), keep=False),
            "SNP",
        ].iloc[0]
        raise ValueError(f"Selected {label} posterior contains duplicate SNP key: {duplicate}")
    if table.loc[:, list(VARIANT_KEYS)].isna().any().any():
        raise ValueError(f"Selected {label} posterior contains missing variant keys")
    for column in ("BETA", "VAR_BETA"):
        values = pd.to_numeric(table[column], errors="coerce").to_numpy(dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError(f"Selected {label} posterior column {column} must be finite")
        if column == "VAR_BETA" and np.any(values <= 0.0):
            raise ValueError(
                f"Selected {label} posterior VAR_BETA must be finite and greater than zero"
            )
        table[column] = values
    return table


def _one_dimensional(values: np.ndarray, label: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional vector")
    return array
