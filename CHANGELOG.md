# Changelog

## 1.0.2 - unreleased

### Corrected

- Defined Stage-1 `BETA` as `gamma*mu` and `VAR_BETA` as the marginal posterior
  variance.
- Kept `var_prior` as the first-E-step variance initialization and restored the
  subsequent EM update of `tau_beta`.
- Added quantitative R2 and binary AUC selection for the fixed 10 by 10 M1/M2
  candidate grid.
- Replaced the reconstructed common-global-mean M3 with directional, pairwise
  base-to-target calibration under `lambda_j ~ Uniform(0,1)`.
- Changed M3 to read one selected target and one selected base posterior,
  consume `VAR_BETA` directly, and emit one calibrated effect vector plus
  convergence diagnostics.
- Replaced the former multi-score ensemble with exactly two target predictors:
  selected target Stage-1 and calibrated Stage-2 scores.
- Added explicit, disjoint target validation and test inputs. Test outcomes are
  read only after Stage-3 fitting and prediction.
- Corrected SuperLearner libraries to Lasso/ridge/neural network for
  quantitative traits and Lasso/neural network with `method.AUC` for binary
  traits.
- Added scientific model identifiers to checkpoint/configuration hashes so old
  M3 and ensemble checkpoints cannot be reused.
- Excluded PLINK2's `NAMED_ALLELE_DOSAGE_SUM` diagnostic from score predictors.
- Added explicit CHR/POS conflict detection for shared SNP identifiers and
  IID-ordered test phenotype matching.
- Removed non-scientific wall-clock timings from the serialized Stage-3 model
  so clean source and wheel installations produce identical model artifacts.
- Made M3 non-convergence a stage failure while retaining its diagnostic table.
- Rejected chromosome, predictor, or phenotype IID-set mismatches instead of
  silently reducing validation/test cohorts through inner joins.
- Removed the reconstructed quantitative-trait covariate residualization so
  Stage 3 follows the stated two-predictor formula exactly.

### Added

- Direct scientific unit tests for Stage-1 posterior moments, M3 directionality
  and ELBO diagnostics, strict posterior alignment, and Stage-3 data isolation.
- Independent target validation/test synthetic genotype and phenotype files.
- Clearly separated selected posteriors, calibrated posterior, score matrices,
  serialized Stage-3 model, predictions, metrics, hyperparameters, and M3
  diagnostics.

### Validation limitation

The implementation is tested against the supplied mathematical specification
and deterministic synthetic examples. It has not yet been numerically compared
with authoritative historical result-generating manuscript outputs.

No package, container, tag, or release has been published.
