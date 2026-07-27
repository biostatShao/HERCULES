"""HERCULES M3 integration matching the recovered ``vi_bayes_paper`` routine."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True, slots=True)
class M3Result:
    beta_target: np.ndarray
    beta_target_variance: np.ndarray
    mu_global: float
    tau_squared: float
    elbo: float
    converged: bool


def vi_bayes_paper(
    mu_target: np.ndarray,
    sigma_target: np.ndarray,
    mu_base: np.ndarray,
    sigma_base: np.ndarray,
    *,
    max_iter: int = 100,
    tol: float = 1e-6,
    sigma0_squared: float = 1e6,
    a: float = 0.001,
    b: float = 0.001,
) -> M3Result:
    """Reproduce the recovered Rcpp M3 update equations.

    The historical implementation passes columns named ``VAR_BETA`` as the
    ``sigma`` arguments and squares them inside this routine.  This behavior is
    intentionally retained until a separate scientific review authorizes a
    variance-interpretation change.
    """

    mu_target = np.asarray(mu_target, dtype=np.float64)
    sigma_target = np.asarray(sigma_target, dtype=np.float64)
    mu_base = np.asarray(mu_base, dtype=np.float64)
    sigma_base = np.asarray(sigma_base, dtype=np.float64)
    if not (mu_target.shape == sigma_target.shape == mu_base.shape == sigma_base.shape):
        raise ValueError("M3 posterior vectors must have identical shapes")
    if mu_target.ndim != 1 or mu_target.size == 0:
        raise ValueError("M3 posterior vectors must be non-empty one-dimensional arrays")
    if not all(
        np.isfinite(values).all()
        for values in (mu_target, sigma_target, mu_base, sigma_base)
    ):
        raise ValueError("M3 posterior vectors must contain only finite values")
    if np.any(sigma_target <= 0.0) or np.any(sigma_base <= 0.0):
        raise ValueError(
            "M3 requires strictly positive VAR_BETA inputs under the current HERCULES "
            "variance-column interpretation"
        )

    target_sigma_squared = np.square(sigma_target)
    base_sigma_squared = np.square(sigma_base)
    target_mean = mu_target.copy()
    target_variance = target_sigma_squared.copy()
    base_mean = mu_base.copy()
    base_variance = base_sigma_squared.copy()

    n_variants = mu_target.size
    mu_global = 0.5 * (float(np.mean(target_mean)) + float(np.mean(base_mean)))
    sigma_global_squared = 1e6
    tau_squared = 1.0
    previous_elbo = -np.inf
    converged = False

    for iteration in range(max_iter):
        target_precision = 1.0 / target_sigma_squared
        tau_precision = 1.0 / tau_squared
        target_mean = (
            target_precision * mu_target + tau_precision * mu_global
        ) / (target_precision + tau_precision)
        target_variance = 1.0 / (target_precision + tau_precision)

        base_precision = 1.0 / base_sigma_squared
        base_mean = (
            base_precision * mu_base + tau_precision * mu_global
        ) / (base_precision + tau_precision)
        base_variance = 1.0 / (base_precision + tau_precision)

        mu_global = (float(np.sum(target_mean)) + float(np.sum(base_mean))) / (
            2.0 * n_variants + 1e-10
        )
        sigma_global_squared = 1.0 / (
            1.0 / sigma0_squared + 2.0 * n_variants / tau_squared
        )
        numerator = (
            float(np.sum(target_variance + np.square(target_mean - mu_global)))
            + float(np.sum(base_variance + np.square(base_mean - mu_global)))
            + 2.0 * b
        )
        tau_squared = numerator / (2.0 * n_variants + 2.0 * a + 2.0)

        tau_precision = 1.0 / tau_squared
        elbo = float(
            np.sum(
                -0.5 * np.log(2.0 * np.pi * target_sigma_squared)
                - 0.5
                * (np.square(mu_target - target_mean) + target_variance)
                / target_sigma_squared
            )
            + np.sum(
                -0.5 * np.log(2.0 * np.pi * base_sigma_squared)
                - 0.5
                * (np.square(mu_base - base_mean) + base_variance)
                / base_sigma_squared
            )
            + np.sum(
                -0.5 * np.log(2.0 * np.pi * tau_squared)
                - 0.5
                * tau_precision
                * (np.square(target_mean - mu_global) + target_variance)
            )
            + np.sum(
                -0.5 * np.log(2.0 * np.pi * tau_squared)
                - 0.5
                * tau_precision
                * (np.square(base_mean - mu_global) + base_variance)
            )
            - 0.5 * np.log(2.0 * np.pi * sigma0_squared)
            - 0.5 * (mu_global**2 + sigma_global_squared) / sigma0_squared
            + a * np.log(b)
            - math.lgamma(a)
            - (a + 1.0) * np.log(tau_squared)
            - b * tau_precision
            + np.sum(0.5 * np.log(2.0 * np.pi * target_variance) + 0.5)
            + np.sum(0.5 * np.log(2.0 * np.pi * base_variance) + 0.5)
        )

        if iteration > 0 and abs(elbo - previous_elbo) < tol:
            converged = True
            break
        previous_elbo = elbo

    return M3Result(
        beta_target=target_mean,
        beta_target_variance=target_variance,
        mu_global=mu_global,
        tau_squared=tau_squared,
        elbo=previous_elbo,
        converged=converged,
    )


def integrate_posterior_tables(
    m1_path: str | Path,
    m2_path: str | Path,
    output_path: str | Path,
    *,
    max_iter: int = 1000,
    tol: float = 1e-6,
) -> Path:
    """Inner-join M1/M2 candidates and write the 100-column M3 score table."""

    keys = ["CHR", "SNP", "POS", "A1", "A2"]
    posterior_columns = [
        *keys,
        *(f"BETA_{index}" for index in range(100)),
        *(f"VAR_BETA_{index}" for index in range(100)),
    ]
    m1 = pd.read_csv(m1_path, sep="\t", usecols=posterior_columns)
    m2 = pd.read_csv(m2_path, sep="\t", usecols=posterior_columns)
    merged = m1.merge(m2, on=keys, how="inner", suffixes=(".x", ".y"), sort=True)
    if merged.empty:
        raise ValueError(f"M1/M2 posterior intersection is empty: {m1_path}, {m2_path}")

    model_indices = sorted(
        int(column.removeprefix("BETA_").removesuffix(".x"))
        for column in merged.columns
        if column.startswith("BETA_") and column.endswith(".x")
    )
    if model_indices != list(range(100)):
        raise ValueError(
            "M3 requires the recovered 100-candidate grid with BETA_0..BETA_99; "
            f"found {len(model_indices)} candidates"
        )

    beta_columns: dict[str, np.ndarray] = {}
    for model_index in model_indices:
        result = vi_bayes_paper(
            merged[f"BETA_{model_index}.x"].to_numpy(),
            merged[f"VAR_BETA_{model_index}.x"].to_numpy(),
            merged[f"BETA_{model_index}.y"].to_numpy(),
            merged[f"VAR_BETA_{model_index}.y"].to_numpy(),
            max_iter=max_iter,
            tol=tol,
        )
        if not np.isfinite(result.beta_target).all():
            raise ValueError(f"M3 candidate {model_index} produced non-finite effects")
        beta_columns[f"BETA_{model_index}"] = result.beta_target

    output = pd.concat(
        [merged.loc[:, keys].reset_index(drop=True), pd.DataFrame(beta_columns)],
        axis=1,
    )

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(destination, sep="\t", index=False)
    return destination
