"""PLINK2 scoring and chromosome aggregation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .process import run_process


def run_plink2_score(
    *,
    plink2: str,
    genotype_prefix: str | Path,
    score_file: str | Path,
    output_prefix: str | Path,
    grid: bool,
) -> Path:
    output = Path(output_prefix)
    output.parent.mkdir(parents=True, exist_ok=True)
    command: list[str | Path] = [
        plink2,
        "--bfile",
        genotype_prefix,
        "--score",
        score_file,
        "2",
        "4",
    ]
    if not grid:
        command.append("6")
    command.extend(["cols=+scoresums,-scoreavgs"])
    if grid:
        command.extend(["--score-col-nums", "6-105"])
    command.extend(["--rm-dup", "force-first", "--out", output])
    result = run_process(command)
    Path(f"{output}.stdout.log").write_text(result.stdout, encoding="utf-8")
    Path(f"{output}.stderr.log").write_text(result.stderr, encoding="utf-8")
    score_path = Path(f"{output}.sscore")
    if not score_path.is_file():
        raise FileNotFoundError(f"PLINK2 completed without creating {score_path}")
    return score_path


def aggregate_score_files(paths: list[Path], output_path: str | Path) -> Path:
    """Inner-join chromosome score files by IID and sum matching score columns."""

    if not paths:
        raise ValueError("At least one PLINK2 score file is required")
    aggregate: pd.DataFrame | None = None
    score_columns: list[str] = []
    for path in paths:
        table = pd.read_csv(path, sep=r"\s+")
        if "IID" not in table.columns:
            raise ValueError(f"PLINK2 score file has no IID column: {path}")
        if table["IID"].duplicated().any():
            duplicate = table.loc[table["IID"].duplicated(keep=False), "IID"].iloc[0]
            raise ValueError(f"PLINK2 score file contains duplicate IID {duplicate}: {path}")
        # PLINK2 may also emit NAMED_ALLELE_DOSAGE_SUM. That diagnostic is not
        # a genetic score and must never enter Stage 3.
        current_scores = [
            column
            for column in table.columns
            if column.startswith("SCORE") and column.endswith("_SUM")
        ]
        if not current_scores:
            raise ValueError(f"PLINK2 score file has no *_SUM columns: {path}")
        current = table.loc[:, ["IID", *current_scores]].copy()
        if aggregate is None:
            aggregate = current
            score_columns = current_scores
            continue
        if current_scores != score_columns:
            raise ValueError(f"PLINK2 score schemas differ across chromosomes: {path}")
        aggregate_iids = set(aggregate["IID"])
        current_iids = set(current["IID"])
        if aggregate_iids != current_iids:
            missing = sorted(str(value) for value in aggregate_iids - current_iids)
            extra = sorted(str(value) for value in current_iids - aggregate_iids)
            raise ValueError(
                "PLINK2 score sample sets differ across chromosomes for "
                f"{path}; missing={len(missing)}"
                f"{f' (example {missing[0]})' if missing else ''}, "
                f"extra={len(extra)}{f' (example {extra[0]})' if extra else ''}"
            )
        joined = aggregate.merge(
            current,
            on="IID",
            how="inner",
            suffixes=(".x", ".y"),
            sort=False,
            validate="one_to_one",
        )
        aggregate = joined.loc[:, ["IID"]].copy()
        for column in score_columns:
            aggregate[column] = joined[f"{column}.x"] + joined[f"{column}.y"]

    assert aggregate is not None
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    aggregate.to_csv(destination, sep="\t", index=False)
    return destination
