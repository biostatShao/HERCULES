"""Dependency diagnostics used by ``hercules doctor``."""

from __future__ import annotations

import importlib
import os
import platform
import shutil
import sys
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from .config import ToolConfig
from .exceptions import ProcessExecutionError
from .process import run_process


@dataclass(frozen=True, slots=True)
class Diagnostic:
    category: str
    name: str
    available: bool
    detail: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def collect_diagnostics(tool_config: ToolConfig | None = None) -> tuple[Diagnostic, ...]:
    configured = tool_config or ToolConfig()
    checks: list[Diagnostic] = [
        Diagnostic(
            "mandatory runtime",
            "Python",
            (3, 10) <= sys.version_info < (3, 13),
            sys.version.split()[0],
        ),
        _module("mandatory runtime", "yaml", "PyYAML"),
        _module("scientific runtime", "numpy", "NumPy"),
        _module("scientific runtime", "scipy", "SciPy"),
        _module("scientific runtime", "pandas", "pandas"),
        _module("scientific runtime", "magenpy", "magenpy"),
        _module("scientific runtime", "hercules.core.model.vi.e_step", "native Cython E-step"),
        _module("scientific runtime", "hercules.core.model.vi.e_step_cpp", "native C++ E-step"),
        _executable("full pipeline", "Rscript", configured.rscript),
        _executable("full pipeline", "PLINK 1.9", configured.plink),
        _executable("full pipeline", "PLINK2", configured.plink2),
        _module("build only", "Cython", "Cython"),
        _executable("build only", "cc", alternatives=("gcc", "clang", "cl")),
        _executable("optional performance", "pkg-config", alternatives=("pkgconf",)),
    ]
    checks.append(_r_packages(configured.rscript))
    checks.append(
        Diagnostic("platform", "system", True, f"{platform.system()} {platform.release()}")
    )
    return tuple(checks)


def mandatory_runtime_ready(checks: Iterable[Diagnostic]) -> bool:
    return all(item.available for item in checks if item.category == "mandatory runtime")


def _module(category: str, module: str, label: str) -> Diagnostic:
    try:
        imported = importlib.import_module(module)
    except (ImportError, ModuleNotFoundError, OSError, ValueError) as exc:
        return Diagnostic(category, label, False, f"not importable: {exc}")
    version = getattr(imported, "__version__", None)
    return Diagnostic(category, label, True, str(version) if version else "available")


def _executable(
    category: str,
    name: str,
    command: str | None = None,
    *,
    alternatives: tuple[str, ...] = (),
) -> Diagnostic:
    candidates = tuple(candidate for candidate in (command or name, *alternatives) if candidate)
    found = None
    for candidate in candidates:
        found = _resolve_executable(candidate)
        if found is not None:
            break
    detail = found or f"configured as {command or name!r}; not found"
    return Diagnostic(category, name, found is not None, detail)


def _resolve_executable(command: str) -> str | None:
    candidate = Path(command)
    if candidate.is_absolute() or candidate.parent != Path("."):
        return str(candidate) if candidate.is_file() and os.access(candidate, os.X_OK) else None
    return shutil.which(command)


def _r_packages(rscript: str) -> Diagnostic:
    resolved = _resolve_executable(rscript)
    packages = ("SuperLearner", "glmnet", "nnet", "pROC", "data.table")
    if resolved is None:
        return Diagnostic("full pipeline", "R packages", False, "Rscript is unavailable")
    package_vector = ",".join(f'"{package}"' for package in packages)
    expression = (
        f"p<-c({package_vector});m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];"
        "if(length(m)){cat(paste(m,collapse=','));quit(status=1)}"
    )
    try:
        run_process([resolved, "-e", expression], timeout=30)
    except ProcessExecutionError as exc:
        return Diagnostic("full pipeline", "R packages", False, str(exc).splitlines()[-1])
    return Diagnostic("full pipeline", "R packages", True, ", ".join(packages))
