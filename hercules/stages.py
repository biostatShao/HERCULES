"""Canonical HERCULES stage definitions."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class StageSpec:
    """The single source of truth for a HERCULES workflow stage."""

    stage_id: str
    display_name: str
    output_prefix: str
    checkpoint_name: str
    dependencies: tuple[str, ...]
    input_requirements: tuple[str, ...]
    output_schema: tuple[str, ...]


STAGES: dict[str, StageSpec] = {
    "m1": StageSpec(
        stage_id="m1",
        display_name="HERCULES M1",
        output_prefix="HERCULES_M1",
        checkpoint_name="m1.complete.json",
        dependencies=(),
        input_requirements=(
            "target summary statistics",
            "target-ancestry LD",
            "target validation genotype and phenotype",
        ),
        output_schema=("CHR", "SNP", "POS", "A1", "A2", "BETA", "VAR_BETA"),
    ),
    "m2": StageSpec(
        stage_id="m2",
        display_name="HERCULES M2",
        output_prefix="HERCULES_M2",
        checkpoint_name="m2.complete.json",
        dependencies=(),
        input_requirements=(
            "base summary statistics",
            "base-ancestry LD",
            "base validation genotype and phenotype",
        ),
        output_schema=("CHR", "SNP", "POS", "A1", "A2", "BETA", "VAR_BETA"),
    ),
    "m3": StageSpec(
        stage_id="m3",
        display_name="HERCULES M3",
        output_prefix="HERCULES_M3",
        checkpoint_name="m3.complete.json",
        dependencies=("m1", "m2"),
        input_requirements=(
            "selected target Stage-1 posterior",
            "selected base Stage-1 posterior",
        ),
        output_schema=(
            "CHR",
            "SNP",
            "POS",
            "A1",
            "A2",
            "BETA",
            "VAR_BETA",
        ),
    ),
    "ensemble": StageSpec(
        stage_id="ensemble",
        display_name="HERCULES ensemble",
        output_prefix="HERCULES_ensemble",
        checkpoint_name="ensemble.complete.json",
        dependencies=("m1", "m3"),
        input_requirements=(
            "target validation genotype and phenotype",
            "target test genotype",
            "selected target Stage-1 posterior",
            "calibrated Stage-2 posterior",
        ),
        output_schema=("IID", "target_stage1_score", "calibrated_stage2_score"),
    ),
}

_ALIASES: dict[str, str] = {
    "m1": "m1",
    "hercules_m1": "m1",
    "m2": "m2",
    "hercules_m2": "m2",
    "m3": "m3",
    "hercules_m3": "m3",
    "ensemble": "ensemble",
    "hercules_ensemble": "ensemble",
}


def get_stage(name: str) -> StageSpec:
    """Resolve a canonical stage name."""

    normalized = name.strip().lower()
    try:
        return STAGES[_ALIASES[normalized]]
    except KeyError as exc:
        valid = ", ".join(STAGES)
        raise KeyError(f"Unknown HERCULES stage {name!r}; expected one of: {valid}") from exc


def execution_order(target: str = "ensemble") -> tuple[StageSpec, ...]:
    """Return the dependency-safe stage order through *target*."""

    target_spec = get_stage(target)
    ordered: list[StageSpec] = []
    seen: set[str] = set()

    def visit(stage: StageSpec) -> None:
        if stage.stage_id in seen:
            return
        for dependency in stage.dependencies:
            visit(STAGES[dependency])
        seen.add(stage.stage_id)
        ordered.append(stage)

    visit(target_spec)
    return tuple(ordered)
