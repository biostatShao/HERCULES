from __future__ import annotations

import sys

from hercules.config import ToolConfig
from hercules.diagnostics import collect_diagnostics


def test_diagnostics_honor_configured_tool_paths() -> None:
    checks = collect_diagnostics(
        ToolConfig(plink=sys.executable, plink2=sys.executable, rscript=sys.executable)
    )
    by_name = {item.name: item for item in checks}
    assert by_name["PLINK 1.9"].available
    assert by_name["PLINK2"].available
    assert by_name["Rscript"].available
