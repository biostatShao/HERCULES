"""Strict FastGWA input validation and internal prior preparation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


FASTGWA_REQUIRED_COLUMNS = (
    "CHR",
    "SNP",
    "POS",
    "A1",
    "A2",
    "N",
    "AF1",
    "BETA",
    "SE",
    "P",
)
VARIANCE_PRIOR_COLUMN = "var_prior"


@dataclass(frozen=True, slots=True)
class PreparedSumstats:
    path: Path
    variants: int
    used_var_prior: bool


def validate_fastgwa_header(path: str | Path) -> tuple[str, ...]:
    """Validate the exact required FastGWA column names without reading rows."""

    input_path = Path(path)
    try:
        columns = tuple(
            pd.read_csv(input_path, sep="\t", compression="infer", nrows=0).columns
        )
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise ValueError(f"Could not read FastGWA header from {input_path}: {exc}") from exc
    missing = [column for column in FASTGWA_REQUIRED_COLUMNS if column not in columns]
    if missing:
        raise ValueError(
            f"FastGWA file {input_path} is missing required columns: "
            + ", ".join(missing)
        )
    return columns


def prepare_fastgwa_sumstats(
    input_path: str | Path,
    output_path: str | Path,
) -> PreparedSumstats:
    """Create the internal FastGWA table used by M1 or M2.

    The user's ``P`` column retains its normal association P-value meaning.
    When ``var_prior`` is present, its positive values become the internal
    per-SNP prior variances. When it is absent, a neutral variance of one is
    used for every SNP. The source file is never modified.
    """

    source = Path(input_path)
    destination = Path(output_path)
    validate_fastgwa_header(source)
    try:
        table = pd.read_csv(source, sep="\t", compression="infer")
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise ValueError(f"Could not read FastGWA file {source}: {exc}") from exc
    if table.empty:
        raise ValueError(f"FastGWA file {source} contains no variants")

    _validate_identifiers(table, source)
    _validate_numeric_columns(table, source)

    used_var_prior = VARIANCE_PRIOR_COLUMN in table.columns
    if used_var_prior:
        prior = _numeric_array(table[VARIANCE_PRIOR_COLUMN], VARIANCE_PRIOR_COLUMN, source)
        invalid = ~np.isfinite(prior) | (prior <= 0.0)
        if np.any(invalid):
            row = int(np.flatnonzero(invalid)[0]) + 2
            raise ValueError(
                f"FastGWA column {VARIANCE_PRIOR_COLUMN!r} in {source} must contain "
                f"finite values greater than zero; invalid value at file row {row}."
            )
    else:
        prior = np.ones(len(table), dtype=np.float64)

    internal = table.loc[:, FASTGWA_REQUIRED_COLUMNS].copy()
    internal.loc[:, "P"] = prior
    destination.parent.mkdir(parents=True, exist_ok=True)
    internal.to_csv(destination, sep="\t", index=False)
    return PreparedSumstats(destination, len(internal), used_var_prior)


def _validate_identifiers(table: pd.DataFrame, source: Path) -> None:
    for column in ("CHR", "SNP", "A1", "A2"):
        if table[column].isna().any() or table[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"FastGWA column {column!r} in {source} contains missing values")


def _validate_numeric_columns(table: pd.DataFrame, source: Path) -> None:
    arrays = {
        column: _numeric_array(table[column], column, source)
        for column in ("POS", "N", "AF1", "BETA", "SE", "P")
    }
    for column, values in arrays.items():
        if not np.isfinite(values).all():
            raise ValueError(f"FastGWA column {column!r} in {source} must contain finite values")
    if np.any(arrays["POS"] <= 0.0):
        raise ValueError(f"FastGWA column 'POS' in {source} must be greater than zero")
    if np.any(arrays["N"] <= 0.0):
        raise ValueError(f"FastGWA column 'N' in {source} must be greater than zero")
    if np.any((arrays["AF1"] < 0.0) | (arrays["AF1"] > 1.0)):
        raise ValueError(f"FastGWA column 'AF1' in {source} must be between zero and one")
    if np.any(arrays["SE"] <= 0.0):
        raise ValueError(f"FastGWA column 'SE' in {source} must be greater than zero")
    if np.any((arrays["P"] < 0.0) | (arrays["P"] > 1.0)):
        raise ValueError(f"FastGWA column 'P' in {source} must be between zero and one")


def _numeric_array(series: pd.Series, column: str, source: Path) -> np.ndarray:
    converted = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
    if np.isnan(converted).any() and not series.isna().any():
        raise ValueError(f"FastGWA column {column!r} in {source} contains non-numeric values")
    return converted
