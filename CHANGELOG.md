# Changelog

## 1.0.1 - unreleased

### Added

- Unified `hercules` Python namespace and command-line interface.
- YAML configuration, validation, stage registry, manifests, and checkpoints.
- Self-contained Python, Cython, and C++ scientific runtime.
- M1, M2, M3, PLINK2 scoring, and R ensemble execution.
- Deterministic quantitative and binary example data.
- Source installation through `pip install .`.
- Linux build, clean-install, and test workflow.
- Explicit stage-specific fixed per-SNP variance priors for M1 and M2. FastGWA
  `P` is interpreted as the precomputed variance and converted to
  `tau_beta = 1 / PVAL` without subsequent M-step replacement.

### Validation

- Native source build and clean installation passed on Linux CPython 3.11.
- All 57 automated tests passed.
- Quantitative and binary end-to-end example workflows passed.
- The corrected editable and clean-wheel quantitative outputs matched exactly
  across M1, M2, M3, scoring, and ensemble tables (maximum absolute numerical
  difference `0.0`).
- The deterministic pre-unification and clean-package example outputs matched
  exactly across M1, M2, M3, scoring, and ensemble tables (maximum absolute
  numerical difference `0.0`); this comparison predates the fixed-prior
  correction and is not evidence of equivalence for the corrected inference.
- One historical M3 comparison passed the documented float64 tolerance.
- Historical M1/M2 comparisons did not pass the predefined float32 tolerance.
- Full scientific equivalence has not yet been demonstrated.

No package, container, tag, or release has been published.
