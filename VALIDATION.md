# Validation status

This document records executed checks only. It does not claim equivalence to
historical manuscript outputs unless such a comparison is explicitly listed.

## Scientific acceptance criteria

- Stage-1 output uses `BETA = gamma * mu` and marginal posterior `VAR_BETA`.
- `var_prior` affects the first E-step only; the M-step then updates
  `tau_beta = pi*M/sum(zeta)`.
- `pi` and residual variance remain fixed within every one of the 100 grid
  candidates.
- M3 is directional, pairwise, has `lambda ~ Uniform(0,1)`, reads selected
  posterior tables, and uses `VAR_BETA` without resquaring.
- Stage 3 contains exactly target Stage-1 and calibrated Stage-2 scores.
- Stage 3 is trained on explicit target validation data and applied unchanged
  to disjoint target test data.
- Quantitative SuperLearner uses Lasso/ridge/neural network; binary uses
  Lasso/neural network with `method.AUC`.

## Automated tests

The suite directly covers the two Stage-1 posterior formulas, initialization
and EM replacement, fixed grid hyperparameters, grid cardinality, directional
M3 behavior, bounded lambda, variance handling, selected-posterior schemas,
alignment failures, finite convergence diagnostics, two-column Stage-3 input,
sample-overlap rejection, and the required learner libraries.

Executed command:

```bash
pytest -q
```

Result: **89 passed**.

## Linux end-to-end smoke tests

Environment:

- Linux x86-64;
- CPython 3.11;
- seed 7209, one inference thread, one parallel job;
- PLINK 1.9, PLINK2;
- R 4.1.2 with SuperLearner, glmnet, nnet, pROC, and data.table.

Both deterministic quantitative and binary fixtures execute the real native
M1/M2 inference, directional M3 integration, PLINK2 validation/test scoring,
and R SuperLearner workflow. Synthetic metrics are smoke-test diagnostics only
and do not estimate performance on real cohorts.

Observed final smoke metrics:

| Fixture | Metric | Value |
|---|---:|---:|
| Quantitative | R2 | 0.984217268991223 |
| Binary | AUC | 0.99937343358396 |

## Installation reproducibility

The 1.0.2 sdist and Linux CPython 3.11 wheel were built in isolation. Each was
installed into a separate environment outside the checkout with its declared
runtime dependencies. Both environments completed the quantitative and binary
workflows under seed 7209.

The following outputs matched byte-for-byte after decompressing gzip tables:

- selected M1 and M2 posteriors;
- selected hyperparameters and validation metrics;
- calibrated M3 posterior and convergence diagnostics;
- validation/test p1 and p2 score matrices;
- Stage-3 metadata, predictions, final metrics, and serialized model object.

The maximum observed absolute and relative numerical differences were both
`0.0`.

## Remaining scientific validation limitation

The previous reconstructed common-global-mean M3 is scientifically incorrect
and is not an acceptance reference. Authoritative historical outputs from the
manuscript's original directional lambda-based Stage 2 have not yet been
provided. Consequently, the current implementation can be tested against the
stated model and independent toy calculations, but historical manuscript-result
equivalence has not yet been demonstrated.
