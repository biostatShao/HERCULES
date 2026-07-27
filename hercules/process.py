"""Safe external-process execution utilities."""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from .exceptions import ProcessExecutionError


@dataclass(frozen=True, slots=True)
class ProcessResult:
    args: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


def run_process(
    args: Sequence[str | os.PathLike[str]],
    *,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
) -> ProcessResult:
    """Run an argument vector, capture logs, and raise an actionable error."""

    argv = tuple(os.fspath(arg) for arg in args)
    if not argv:
        raise ValueError("External command cannot be empty")

    try:
        completed = subprocess.run(
            argv,
            cwd=Path(cwd) if cwd is not None else None,
            env=dict(env) if env is not None else None,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
            shell=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ProcessExecutionError(
            f"Command timed out after {timeout} seconds: {argv!r}"
        ) from exc
    except OSError as exc:
        location = f" in {Path(cwd)}" if cwd is not None else ""
        raise ProcessExecutionError(f"Could not start command{location}: {argv!r}: {exc}") from exc
    result = ProcessResult(argv, completed.returncode, completed.stdout, completed.stderr)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "no process output"
        raise ProcessExecutionError(
            f"Command failed with exit code {completed.returncode}: {argv!r}\n{detail}"
        )
    return result
