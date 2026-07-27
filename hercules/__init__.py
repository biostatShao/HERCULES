"""Public Python interface for the packaged HERCULES workflow.

The top-level import remains lightweight; scientific classes and native
extensions are loaded only when their implementation modules are imported.
"""

from .stages import STAGES, StageSpec, get_stage

__all__ = ["STAGES", "StageSpec", "__version__", "get_stage"]

__version__ = "0.1.0.dev0"
