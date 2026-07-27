"""Safety gates separating baseline recovery from public scientific execution."""

from __future__ import annotations

import os

from .exceptions import BaselineUnavailableError

_REFERENCE_TOKEN = "I_ACKNOWLEDGE_UNVALIDATED_REFERENCE_EXECUTION"


def require_reference_execution_authorized() -> None:
    """Block inherited model construction outside an explicit recovery harness."""

    if os.environ.get("HERCULES_REFERENCE_HARNESS") == _REFERENCE_TOKEN:
        return
    raise BaselineUnavailableError(
        "This reference-only inference path is baseline-gated. It may only be loaded by an "
        "explicit recovery harness before differential fixtures are frozen."
    )
