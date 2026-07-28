# Included example data

This directory contains a deterministic, public smoke-test fixture:

- two synthetic ancestries;
- 160 samples per ancestry;
- 64 variants on chromosome 22;
- PLINK genotype files;
- windowed LD stores;
- summary statistics;
- separate target validation/test genotypes and phenotypes;
- ancestry-matched quantitative and binary Stage-1 validation phenotypes;
- functional annotations and per-SNP values;
- ready-to-run quantitative and binary YAML configurations.

In both FastGWA summary-statistics files, `P` is a normal association P value
and `var_prior` is a positive synthetic per-SNP initial variance. M1 uses the
target file and M2 uses the base file independently.

HERCULES converts these values to the initial precision
`tau_beta_j = 1/var_prior_j`. The EM M-step subsequently updates `tau_beta`;
the input values are not fixed throughout fitting.

M3 reads only the selected target/base posterior means and marginal variances,
then produces one directional calibrated target effect vector. The final
learner is fitted on the target validation cohort with exactly the target
Stage-1 and calibrated Stage-2 scores and applied to the disjoint target test
cohort.

The fixture contains no private or real participant data. Run commands from the
repository root so that the relative paths in the YAML files resolve correctly.

```bash
hercules run --config examples/data/hercules.quantitative.yaml
hercules run --config examples/data/hercules.binary.yaml
```

Generated outputs are written under `examples/output/` and are ignored by Git.
