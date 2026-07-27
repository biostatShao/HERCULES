from __future__ import annotations

from pathlib import Path

from hercules.resources import ensemble_script_path, example_config_path, r_interface_path


def test_source_resources_are_locatable() -> None:
    assert r_interface_path().name == "HERCULES.R"
    assert example_config_path().name == "hercules.example.yaml"
    assert ensemble_script_path().name == "ensemble_superlearner.R"


def test_packaged_resources_match_checkout_interfaces() -> None:
    root = Path(__file__).resolve().parent.parent
    assert r_interface_path().read_text(encoding="utf-8").rstrip() == (
        root / "R" / "HERCULES.R"
    ).read_text(encoding="utf-8").rstrip()
    assert example_config_path().read_text(encoding="utf-8").rstrip() == (
        root / "examples" / "hercules.example.yaml"
    ).read_text(encoding="utf-8").rstrip()
