from __future__ import annotations

import pytest

from hercules.stages import STAGES, execution_order, get_stage


@pytest.mark.parametrize(
    ("alias", "stage_id"),
    [
        ("HERCULES_M1", "m1"),
        ("HERCULES_M2", "m2"),
        ("HERCULES_M3", "m3"),
        ("HERCULES_ensemble", "ensemble"),
    ],
)
def test_canonical_stage_aliases(alias: str, stage_id: str) -> None:
    assert get_stage(alias).stage_id == stage_id


def test_stage_registry_has_canonical_names() -> None:
    assert [STAGES[key].output_prefix for key in ("m1", "m2", "m3", "ensemble")] == [
        "HERCULES_M1",
        "HERCULES_M2",
        "HERCULES_M3",
        "HERCULES_ensemble",
    ]


def test_m3_dependency_order() -> None:
    assert [stage.stage_id for stage in execution_order("m3")] == ["m1", "m2", "m3"]
