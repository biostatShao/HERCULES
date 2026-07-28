# HERCULES

HERCULES is a cross-ancestry polygenic risk score workflow with three model
stages and a final ensemble:

- M1 fits the target-ancestry model grid and scores the target genotype.
- M2 fits the base-ancestry model grid and scores the target genotype.
- M3 integrates the paired M1 and M2 posterior tables.
- The ensemble combines selected and grid scores for quantitative or binary
  traits.

All public interfaces use one name:

```text
Python package: hercules
Command line:   hercules
R function:     HERCULES()
Configuration:  hercules.yaml
Outputs:        HERCULES_M1, HERCULES_M2, HERCULES_M3, HERCULES_ensemble
```

## Supported platform

The validated platform is Linux x86-64 with Python 3.11. HERCULES contains
Cython/C++/OpenMP extensions, so installation from source requires a C/C++
compiler. Native Windows installation has not been validated; Windows users
should use WSL2, a Linux server, or a Linux container.

The complete workflow also needs:

- PLINK 1.9;
- PLINK2;
- R with `data.table`, `SuperLearner`, `glmnet`, and `pROC`;
- legal user-provided GWAS, LD, genotype, phenotype, and annotation inputs.

## Installation

### Download or clone, then install

This is the recommended installation method because it compiles the native
extensions for the user's own Linux environment instead of relying on a wheel
built on another computer.

```bash
git clone https://github.com/biostatShao/HERCULES.git
cd HERCULES

python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

Development installation:

```bash
python -m pip install -e ".[test]"
pytest
```

Direct installation from GitHub is also possible after the repository is
published:

```bash
python -m pip install "git+https://github.com/biostatShao/HERCULES.git"
```

Prebuilt wheels may be added later for individually tested Python/platform
combinations. A single wheel is not portable across arbitrary operating
systems, Python versions, or CPU architectures.

## Verify the installation

```bash
hercules --version
hercules --help
hercules doctor
```

With a configured workflow:

```bash
hercules doctor --config hercules.yaml
hercules config validate hercules.yaml
```

## Run the included example

The repository contains a small deterministic two-ancestry dataset under
`examples/data/`. It contains 160 samples per ancestry and 64 variants on
chromosome 22. It contains no private research data.

The example configuration uses paths relative to the repository root. Update
the three executable paths in the YAML files if PLINK, PLINK2, or Rscript are
not on `PATH`.

```bash
hercules config validate examples/data/hercules.quantitative.yaml
hercules run --config examples/data/hercules.quantitative.yaml

hercules config validate examples/data/hercules.binary.yaml
hercules run --config examples/data/hercules.binary.yaml
```

To regenerate the example instead of using the committed fixture:

```bash
python examples/synthetic/generate_fixture.py \
  --output examples/data \
  --plink /path/to/plink \
  --plink2 /path/to/plink2 \
  --rscript /path/to/Rscript
```

## Configure a real analysis

Copy the template and replace its paths:

```bash
cp examples/hercules.example.yaml hercules.yaml
hercules config validate hercules.yaml
hercules run --config hercules.yaml
```

Individual dependency-aware stages are available:

```bash
hercules stage m1 --config hercules.yaml
hercules stage m2 --config hercules.yaml
hercules stage m3 --config hercules.yaml
hercules ensemble --config hercules.yaml
```

Requesting M3 or the ensemble automatically runs missing prerequisite stages.
Checkpoints include the trait, chromosome, stage, and configuration hash.

## Input files

### GWAS summary statistics

M1 reads the target-ancestry summary-statistics file and M2 reads the
base-ancestry file. Both files must be tab-delimited FastGWA-style tables with
the following exact, case-sensitive column names:

| Column | Required | Meaning and validation |
|---|---:|---|
| `CHR` | yes | Chromosome identifier. |
| `SNP` | yes | Variant identifier, normally an rsID. |
| `POS` | yes | Positive base-pair position. |
| `A1` | yes | Effect allele corresponding to `BETA`. |
| `A2` | yes | Other allele. |
| `N` | yes | Positive per-SNP GWAS sample size. |
| `AF1` | yes | Effect-allele frequency in the range 0–1. |
| `BETA` | yes | GWAS effect estimate. |
| `SE` | yes | Positive standard error of `BETA`. |
| `P` | yes | Normal association P value in the range 0–1. |
| `var_prior` | no | Precomputed positive, finite per-SNP effect-size prior variance. |

Example:

```text
CHR  SNP       POS      A1  A2  N      AF1   BETA    SE     P       var_prior
22   rs10001   1000000  G   A   50000  0.15  0.006   0.020  0.7642  0.0200
22   rs10002   1010000  T   C   50000  0.17 -0.007   0.021  0.7389  0.0205
```

`P` is always a statistical P value in the public input. HERCULES does not
modify the user's file. Immediately before M1 or M2 inference it creates an
internal temporary table:

- when `var_prior` is present, its values are used as the fixed per-SNP prior
  variances;
- when `var_prior` is absent, every SNP receives prior variance `1`;
- the model uses prior precision `tau_beta_j = 1 / var_prior_j` and keeps it
  fixed during initialization and every M-step.

If `var_prior` is present, every row must contain a numeric value greater than
zero. Missing, zero, negative, infinite, or non-numeric values cause the run to
stop before model fitting. Additional input columns are permitted but are not
passed to the inference parser.

### LD reference

`inputs.ld_reference.base` and `inputs.ld_reference.target` must point to
compatible magenpy LD stores. A path may contain `{chrom}` or `{chromosome}`;
the placeholder is replaced separately for every configured chromosome. SNP
identifiers, positions and alleles must be compatible with the corresponding
summary-statistics file.

### Genotype files

Genotypes must use a complete PLINK prefix:

- BED format: `.bed`, `.bim` and `.fam`; or
- PGEN format: `.pgen`, `.pvar` and `.psam`.

`inputs.genotype_prefixes.target` is used for final PRS scoring.
`inputs.validation_genotype` is the default M1 validation genotype.
`inputs.genotype_prefixes.base_validation` is the default M2 validation
genotype and may be overridden by `m2.validation_genotype`.

### Phenotype, covariates and validation files

`inputs.phenotype_file` is a tab-delimited table used by the final ensemble. It
must contain `IID`, the column named by `inputs.phenotype_column`, and every
column listed under `inputs.covariates`.

The M1/M2 `validation_phenotype` files use the standard PLINK three-column,
header-free layout:

```text
FID  IID  phenotype
```

An optional `validation_keep` file contains `FID` and `IID`, without a header.
For binary traits, phenotype values must follow the coding expected by the
configured PLINK validation data and the final R ensemble.

### Functional input files

`inputs.functional_annotation` and `inputs.per_snp_heritability` record and
validate the source functional files used to produce `var_prior`. The current
runtime does not estimate `var_prior` from raw annotations; that preprocessing
must be completed before running HERCULES.

## YAML configuration reference

Start from `examples/hercules.example.yaml`. Paths may be absolute or relative
to the directory from which the command is run.

### Analysis identity

| Parameter | Description |
|---|---|
| `trait_name` | Trait identifier used in manifests and checkpoints. It must be safe for use in filenames. |
| `chromosomes` | List of chromosomes to process, for example `[1, 2, 22]`. |
| `base_ancestry` | Label for the base ancestry used by M2. |
| `target_ancestry` | Label for the target ancestry used by M1 and final scoring. |

### `inputs`

| Parameter | Description |
|---|---|
| `summary_statistics.base_path` | Base-ancestry FastGWA file or chromosome path template. |
| `summary_statistics.target_path` | Target-ancestry FastGWA file or chromosome path template. |
| `summary_statistics.base_columns` | Keep `{}` for the strict FastGWA interface. |
| `summary_statistics.target_columns` | Keep `{}` for the strict FastGWA interface. |
| `functional_annotation` | Optional raw annotation file retained as input metadata. Use `""` if unavailable. |
| `per_snp_heritability` | Optional file from which `var_prior` was prepared. Use `""` if unavailable. |
| `ld_reference.base` | Base-ancestry magenpy LD path or template. |
| `ld_reference.target` | Target-ancestry magenpy LD path or template. |
| `genotype_prefixes.base_validation` | Base-ancestry PLINK validation prefix used by M2. |
| `genotype_prefixes.target` | Target-ancestry PLINK prefix used for M1/M2/M3 scoring. |
| `validation_genotype` | Default target validation PLINK prefix used by M1. |
| `phenotype_file` | Final ensemble phenotype/covariate table. |
| `phenotype_column` | Outcome column in `phenotype_file`. |
| `covariates` | Covariate column names; use `[]` when no covariates are required. |
| `trait_type` | `quantitative` or `binary`. |

### Output and executables

| Parameter | Description |
|---|---|
| `output_dir` | Final scores, ensemble results, manifests and checkpoints. |
| `temporary_dir` | Per-chromosome posteriors, internal summary statistics and subprocess files. |
| `tools.plink` | PLINK 1.9 executable name or absolute path. |
| `tools.plink2` | PLINK2 executable name or absolute path. |
| `tools.rscript` | Rscript executable name or absolute path. |

### `execution`

| Parameter | Description |
|---|---|
| `threads` | Native inference threads. Use `1` for the validated deterministic mode. |
| `parallel_jobs` | Number of chromosome/model worker processes. Use `1` for deterministic validation. |
| `seed` | Python and ensemble random seed; the default validated seed is `7209`. |

### `m1` and `m2`

Both stages use the same keys, but M1 reads target-ancestry inputs while M2
reads base-ancestry inputs.

| Parameter | Description |
|---|---|
| `hyperparameter_search` | Use `grid` for the validated workflow. |
| `grid_metric` | Model-selection criterion: normally `validation`; `ELBO` and `pseudo_validation` are supported by the internal runner where their required inputs are available. |
| `pi_steps` | Number of causal-proportion grid values; validated default `10`. |
| `sigma_epsilon_steps` | Number of residual-variance grid values; validated default `10`. |
| `max_iter` | Maximum variational inference iterations; validated default `500`. |
| `sumstats_format` | Must be `fastgwa`. |
| `backend` | Optional magenpy genotype backend; validated value `plink`. |
| `validation_genotype` | Optional stage-specific validation PLINK prefix. M1/M2 have ancestry-specific defaults described above. |
| `validation_phenotype` | PLINK-format phenotype used for grid-model selection. |
| `validation_keep` | Optional PLINK-format sample keep file. |

The fixed prior conversion is automatic and does not require a YAML parameter.

### `m3`, `ensemble`, checkpoints and logging

| Parameter | Description |
|---|---|
| `m3.max_iter` | Maximum M3 integration iterations; validated default `1000`. |
| `m3.tol` | M3 convergence tolerance; validated default `1e-6`. |
| `ensemble.quantitative_learners` | Documents the validated quantitative SuperLearner candidates, `SL.glmnet` and `SL.ridge`. |
| `ensemble.binary_selection` | Documents the validated binary rule, `best_individual_auc`. |
| `checkpoint.enabled` | Write stage completion markers. |
| `checkpoint.resume` | Reuse outputs only when the checkpoint and configuration hash match. |
| `logging.level` | Logging level such as `INFO` or `DEBUG`. |
| `logging.file` | Optional log path; use `""` for console/default logging. |

Supported path placeholders are `{chrom}`, `{chromosome}` and `{ancestry}`.

Configuration precedence is:

```text
CLI override > HERCULES__... environment variable > YAML > package default
```

## R interface

```r
source("R/HERCULES.R")
HERCULES("hercules.yaml")
```

Python owns configuration, orchestration, manifests, checkpoints, and external
process handling. R is used only for the final validated ensemble procedure.

## Scientific validation status

The implementation is executable end to end on the validated Linux platform.
The corrected package passed 66 automated tests. Quantitative and binary
examples completed end to end, and the quantitative example also completed
from a clean wheel installation outside the source checkout. A table-by-table
comparison between the corrected editable installation and corrected wheel
installation covered M1/M2 scores and posteriors, M3, and ensemble outputs;
every compared numeric value had a maximum absolute difference of `0.0`.

The corrected quantitative example completed in 68.94 seconds with 249,556
KiB peak resident memory and produced R2 = 0.0577308741972629. The corrected
binary example completed in 19.32 seconds with 216,928 KiB peak resident memory
and produced AUC = 0.683862433862434. These values validate packaging and
deterministic execution, not scientific performance on real data.

The final reported mean `tau_beta` matched `mean(1/var_prior)` from the stage-specific
input to float32 output precision: absolute differences were approximately
`8.2e-7` for M1 and `1.1e-6` for M2.

One historical M3 comparison passed at `rtol=1e-8, atol=1e-10`. Historical M1
and M2 replays retained exact schemas and SNP order but did not pass the
predefined float32 tolerance. Raw annotation-to-per-SNP preprocessing remains
an upstream input-preparation responsibility. The example FastGWA files use
normal association P values in `P` and positive synthetic variances in
`var_prior`.

The recovered M3 source currently contains no ancestry-bridging `lambda`
parameter and no Beta or Uniform prior for such a parameter. Therefore no
unverified lambda implementation has been added. Aligning a manuscript
`lambda ~ Uniform(0,1)` description with executable code requires the original
analysis source (including its lambda updates) or a separately authorized and
numerically validated scientific implementation.

Therefore:

> The refactor was designed to preserve the identified scientific defaults,
> but full scientific equivalence has not yet been demonstrated.

See [VALIDATION.md](VALIDATION.md) for the executed checks and their limits.

## License and attribution

HERCULES is distributed under the MIT License. Required upstream copyright and
derivative-work attribution are preserved in [LICENSE](LICENSE) and
[NOTICE](NOTICE).
