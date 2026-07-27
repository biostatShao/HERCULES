from __future__ import annotations

import sys

import pytest

from hercules.exceptions import ProcessExecutionError
from hercules.process import run_process


def test_process_uses_argument_vector() -> None:
    result = run_process([sys.executable, "-c", "print('path with spaces')"])
    assert result.returncode == 0
    assert result.stdout.strip() == "path with spaces"


def test_process_failure_includes_exit_code_and_stderr() -> None:
    with pytest.raises(ProcessExecutionError, match="exit code 7") as error:
        run_process([sys.executable, "-c", "import sys; print('failed', file=sys.stderr); sys.exit(7)"])
    assert "failed" in str(error.value)


def test_missing_executable_is_actionable() -> None:
    with pytest.raises(ProcessExecutionError, match="Could not start command"):
        run_process(["definitely-not-a-hercules-executable-93c78"])


def test_timeout_is_actionable() -> None:
    with pytest.raises(ProcessExecutionError, match="timed out"):
        run_process([sys.executable, "-c", "import time; time.sleep(2)"], timeout=0.01)
