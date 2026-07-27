"""Locate installed non-Python interfaces and examples."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path


def r_interface_path() -> Path:
    """Return the installed or source-checkout path to the canonical R interface."""

    return _resource_path("R", "HERCULES.R")


def example_config_path() -> Path:
    """Return the installed or source-checkout example YAML path."""

    return _resource_path("examples", "hercules.example.yaml")


def ensemble_script_path() -> Path:
    """Return the installed R ensemble implementation."""

    return _resource_path("R", "ensemble_superlearner.R")


def _resource_path(directory: str, filename: str) -> Path:
    resource = files("hercules").joinpath("data", directory, filename)
    if resource.is_file():
        return Path(str(resource))
    raise FileNotFoundError(f"HERCULES resource was not installed: {directory}/{filename}")
