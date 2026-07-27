from __future__ import annotations

import os
from pathlib import Path

from hercules.config import HerculesConfig, ToolConfig
from hercules.workflow import _fit_command_prefix, _scientific_subprocess_environment


def test_internal_fit_uses_python_safe_path_mode() -> None:
    assert _fit_command_prefix()[1:] == ("-P", "-m", "hercules.fit_cli")


def test_configured_plink_directories_are_exposed_to_nested_tools(
    tmp_path: Path, monkeypatch
) -> None:
    plink_dir = tmp_path / "plink tools"
    plink2_dir = tmp_path / "plink2 tools"
    plink_dir.mkdir()
    plink2_dir.mkdir()
    plink = plink_dir / "plink"
    plink2 = plink2_dir / "plink2"
    plink.write_bytes(b"")
    plink2.write_bytes(b"")
    monkeypatch.setenv("PATH", "existing-path")
    config = HerculesConfig(tools=ToolConfig(plink=str(plink), plink2=str(plink2)))

    environment = _scientific_subprocess_environment(config)

    path_entries = environment["PATH"].split(os.pathsep)
    assert path_entries[:2] == [str(plink_dir.resolve()), str(plink2_dir.resolve())]
    assert path_entries[-1] == "existing-path"
