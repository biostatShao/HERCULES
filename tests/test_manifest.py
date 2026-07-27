from __future__ import annotations

import json
from pathlib import Path

import pytest

from hercules.manifest import checkpoint_matches, checkpoint_path, configuration_hash
from hercules.stages import get_stage


def test_configuration_hash_is_order_independent() -> None:
    assert configuration_hash({"a": 1, "b": 2}) == configuration_hash({"b": 2, "a": 1})


def test_checkpoint_is_trait_chromosome_stage_and_config_specific(tmp_path: Path) -> None:
    marker = checkpoint_path(
        tmp_path, stage=get_stage("m1"), trait="height", chromosome=1
    )
    marker.parent.mkdir(parents=True)
    marker.write_text(
        json.dumps({"status": "complete", "configuration_hash": "expected"}),
        encoding="utf-8",
    )
    assert "height" in marker.parts
    assert "chr1" in marker.parts
    assert checkpoint_matches(marker, "expected")
    assert not checkpoint_matches(marker, "different")


def test_checkpoint_path_rejects_traversal(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Unsafe trait"):
        checkpoint_path(tmp_path, stage=get_stage("m1"), trait="../height", chromosome=1)


def test_non_mapping_checkpoint_is_not_a_match(tmp_path: Path) -> None:
    marker = tmp_path / "bad.json"
    marker.write_text("[]", encoding="utf-8")
    assert not checkpoint_matches(marker, "expected")
