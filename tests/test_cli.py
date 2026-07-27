from __future__ import annotations

import json
from pathlib import Path

import pytest

from hercules.cli import main


def test_config_validate_without_path_checks(config_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["config", "validate", str(config_path), "--no-check-paths"]) == 0
    assert "Configuration is valid" in capsys.readouterr().out


def test_dry_run_reports_dependency_safe_plan(
    config_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["run", "--config", str(config_path), "--dry-run"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["stages"] == ["m1", "m2", "m3", "ensemble"]


def test_scientific_execution_reports_failed_external_runtime(config_path: Path) -> None:
    with pytest.raises(SystemExit) as error:
        main(["run", "--config", str(config_path)])
    assert error.value.code == 2


def test_doctor_accepts_configured_tool_paths(
    config_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    assert main(["doctor", "--config", str(config_path), "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    names = {item["name"] for item in payload}
    assert {"Rscript", "PLINK 1.9", "PLINK2", "R packages"} <= names
