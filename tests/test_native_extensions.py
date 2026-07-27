from __future__ import annotations

from hercules.diagnostics import collect_diagnostics


def test_native_extension_is_loadable_or_reported_missing() -> None:
    diagnostic = next(
        item for item in collect_diagnostics() if item.name == "native Cython E-step"
    )
    if diagnostic.available:
        from hercules.core.model.vi import e_step

        assert e_step is not None
    else:
        assert "not" in diagnostic.detail.lower()
