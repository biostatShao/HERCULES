# Deterministic synthetic workflow

This generator creates a public, deterministic two-ancestry smoke dataset with
PLINK genotype files, magenpy LD stores, fastGWA-style summary statistics,
functional annotations, per-SNP weights, validation phenotypes, and both
quantitative and binary HERCULES configurations. It contains no private data.
The generated FastGWA files contain normal association P values in `P` and
stage-specific per-SNP prior variances in `var_prior`.

Requirements: Python 3.11 with HERCULES installed, PLINK 1.9, PLINK2, and Rscript.

```bash
python examples/synthetic/generate_fixture.py \
  --output examples/synthetic/data \
  --plink /path/to/plink \
  --plink2 /path/to/plink2 \
  --rscript /path/to/Rscript

hercules config validate examples/synthetic/data/hercules.quantitative.yaml
hercules run --config examples/synthetic/data/hercules.quantitative.yaml
hercules run --config examples/synthetic/data/hercules.binary.yaml
```

The generator fixes the Python RNG seed, creates 160 samples per ancestry and
64 variants on chromosome 22, and configures one thread, one parallel job, and
seed 7209. Generated `data/` and `results/` directories are ignored by Git.
