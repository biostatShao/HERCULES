# Included example data

This directory contains a deterministic, public smoke-test fixture:

- two synthetic ancestries;
- 160 samples per ancestry;
- 64 variants on chromosome 22;
- PLINK genotype files;
- windowed LD stores;
- summary statistics;
- quantitative and binary validation phenotypes;
- functional annotations and per-SNP values;
- ready-to-run quantitative and binary YAML configurations.

In both FastGWA summary-statistics files, the physical `P` column contains a
positive synthetic per-SNP prior variance, not an association P-value. M1 uses
the target file and M2 uses the base file independently.

The fixture contains no private or real participant data. Run commands from the
repository root so that the relative paths in the YAML files resolve correctly.

```bash
hercules run --config examples/data/hercules.quantitative.yaml
hercules run --config examples/data/hercules.binary.yaml
```

Generated outputs are written under `examples/output/` and are ignored by Git.
