# Changelog

## 0.1.0.dev0 - unreleased

### Added

- Unified `hercules` Python namespace and command-line interface.
- YAML configuration, validation, stage registry, manifests, and checkpoints.
- Self-contained Python, Cython, and C++ scientific runtime.
- M1, M2, M3, PLINK2 scoring, and R ensemble execution.
- Deterministic quantitative and binary example data.
- Source installation through `pip install .`.
- Linux build, clean-install, and test workflow.

### Validation

- Native source build and clean installation passed on Linux CPython 3.11.
- All 48 automated tests passed.
- Quantitative and binary end-to-end example workflows passed.
- The deterministic pre-unification and clean-package example outputs matched
  exactly across M1, M2, M3, scoring, and ensemble tables (maximum absolute
  numerical difference `0.0`).
- One historical M3 comparison passed the documented float64 tolerance.
- Historical M1/M2 comparisons did not pass the predefined float32 tolerance.
- Full scientific equivalence has not yet been demonstrated.

No package, container, tag, or release has been published.
