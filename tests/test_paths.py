from __future__ import annotations

from pathlib import Path

import pytest

from hercules.paths import stage_output_path


def test_canonical_output_path_supports_parent_spaces(tmp_path: Path) -> None:
    output = stage_output_path(tmp_path / "results with spaces", trait="height", stage="m1", chromosome=1)
    assert output.name == "HERCULES_M1.chr1.tsv.gz"
    assert output.parent.name == "height"


def test_path_component_traversal_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unsafe trait"):
        stage_output_path("out", trait="../height", stage="m1")
