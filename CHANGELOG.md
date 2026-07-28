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
- Explicit stage-specific per-SNP variance initialization for M1 and M2.
  FastGWA input keeps the normal association P value in `P` and accepts an
  initial variance in optional `var_prior`. Missing `var_prior` values default
  to one. The first E-step uses `tau_beta_j = 1/var_prior_j`; subsequent
  M-steps update `tau_beta` through the original EM equation.

### Validation

- Native source build and clean installation passed on Linux CPython 3.11.
- All 67 automated tests passed, including strict FastGWA input validation,
  both `var_prior` initialization branches, and a direct M-step `tau_beta`
  update test.
- Quantitative and binary end-to-end example workflows passed.
- The corrected editable and clean-wheel quantitative outputs matched exactly
  across M1, M2, M3, scoring, and ensemble tables (maximum absolute numerical
  difference `0.0`).
- A complete quantitative run without `var_prior` verified the documented
  all-ones initialization fallback through M1, M2, M3, scoring, and ensemble
  execution.
- The deterministic pre-unification and clean-package example outputs matched
  exactly across M1, M2, M3, scoring, and ensemble tables (maximum absolute
  numerical difference `0.0`); this comparison predates the fixed-prior
  correction and is not evidence of equivalence for the corrected inference.
- One historical M3 comparison passed the documented float64 tolerance.
- Historical M1/M2 comparisons did not pass the predefined float32 tolerance.
- Full scientific equivalence has not yet been demonstrated.

No package, container, tag, or release has been published.
