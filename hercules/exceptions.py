"""HERCULES-specific exceptions."""


class HerculesError(Exception):
    """Base class for actionable public errors."""


class ConfigurationError(HerculesError):
    """Raised when a configuration is invalid."""


class ProcessExecutionError(HerculesError):
    """Raised when an external process exits unsuccessfully."""


class BaselineUnavailableError(HerculesError):
    """Raised when execution would cross the unresolved numerical baseline gate."""

