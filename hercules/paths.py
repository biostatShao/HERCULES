"""Central output-path generation for HERCULES stages."""

from __future__ import annotations

from pathlib import Path

from .stages import StageSpec, get_stage


def stage_output_path(
    output_dir: str | Path,
    *,
    trait: str,
    stage: str | StageSpec,
    chromosome: str | int | None = None,
    suffix: str = ".tsv.gz",
) -> Path:
    spec = get_stage(stage) if isinstance(stage, str) else stage
    safe_trait = safe_path_component(trait, "trait")
    chromosome_part = (
        "" if chromosome is None else f".chr{safe_path_component(str(chromosome), 'chromosome')}"
    )
    return Path(output_dir) / safe_trait / f"{spec.output_prefix}{chromosome_part}{suffix}"


def safe_path_component(value: str, label: str) -> str:
    invalid_characters = set('<>:"|?*')
    if (
        not value
        or value in {".", ".."}
        or "/" in value
        or "\\" in value
        or any(character in invalid_characters or ord(character) < 32 for character in value)
        or value.endswith((" ", "."))
    ):
        raise ValueError(f"Unsafe {label} path component: {value!r}")
    return value
