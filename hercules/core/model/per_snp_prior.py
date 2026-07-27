"""Validation and conversion for fixed per-SNP effect-size priors."""

from __future__ import annotations

from typing import Any

import numpy as np


def prepare_fixed_per_snp_precision(
    values: Any,
    *,
    input_type: str,
    column: str,
    chromosome: str | int,
    float_precision: str,
    order: str = "F",
) -> np.ndarray:
    """Return an ``(n_snps, 1)`` array of fixed prior precisions.

    HERCULES uses ``tau_beta`` as a precision (inverse variance).  The public
    workflow stores the precomputed per-SNP variance in the FastGWA ``P``
    column, which magenpy exposes to the model as ``PVAL``.
    """

    array = np.asarray(values)
    if array.ndim != 1:
        raise ValueError(
            f"Per-SNP prior column {column!r} for chromosome {chromosome} "
            f"must be one-dimensional; observed shape {array.shape}."
        )
    try:
        array = array.astype(float_precision, copy=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Per-SNP prior column {column!r} for chromosome {chromosome} "
            "must contain numeric values."
        ) from exc

    invalid = ~np.isfinite(array) | (array <= 0.0)
    if np.any(invalid):
        first = int(np.flatnonzero(invalid)[0])
        raise ValueError(
            f"Per-SNP prior column {column!r} for chromosome {chromosome} "
            "must contain finite values greater than zero; "
            f"row {first} contains {array[first]!r}."
        )

    normalized_type = input_type.lower()
    if normalized_type == "variance":
        precision = np.reciprocal(array)
    elif normalized_type == "precision":
        precision = array
    else:
        raise ValueError(
            "Per-SNP prior input_type must be 'variance' or 'precision'; "
            f"received {input_type!r}."
        )

    return np.asarray(precision[:, None], dtype=float_precision, order=order)
