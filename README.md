# HERCULES

HERCULES is a Linux-oriented cross-ancestry polygenic-score workflow. Version
1.0.2 implements three scientific stages:

1. **M1** fits and selects the target-population Stage-1 model.
2. **M2** fits and selects one designated base-population Stage-1 model.
3. **M3** directionally calibrates the selected base posterior to the selected
   target posterior using a SNP-specific ancestry-bridging parameter.

The final Stage-3 learner uses exactly two target-sample scores: the selected
target Stage-1 score and the calibrated Stage-2 score.

The package includes the Python/Cython/C++ inference runtime, workflow CLI,
PLINK2 scoring integration, an R SuperLearner boundary, example configuration,
and deterministic public synthetic data.

## Scientific workflow

### M1 and M2

M1 uses target GWAS statistics, target LD, and target validation samples. M2
uses base GWAS statistics, base LD, and base validation samples. Each stage
fits the existing 10 by 10 grid of fixed `pi` and residual-variance values and
selects one candidate using validation R2 for quantitative traits or validation
AUC for binary traits.

For SNP `j`, optional `var_prior` initializes the first variational iteration:

```text
v_j^(0)        = var_prior_j
tau_beta,j^(0) = 1 / var_prior_j
```

If `var_prior` is absent, the initial variance is one. The values are only
initialization values. After the first E-step, the M-step replaces them using:

```text
zeta_j  = gamma_j * (v_j + mu_j^2)
tau_beta = pi * M / sum_j(zeta_j)
```

`pi` and residual variance stay fixed inside each grid candidate. Stage-1
posterior output is defined as:

```text
BETA_j     = gamma_j * mu_j
VAR_BETA_j = gamma_j * (v_j + mu_j^2) - (gamma_j * mu_j)^2
```

`VAR_BETA` is the marginal posterior variance.

### Directional M3 calibration

For one base/target pair, M3 reads only the two selected Stage-1 posterior
tables. For each aligned SNP:

```text
eta_j | lambda_j, b_base_j, V_base_j
  ~ Normal(lambda_j * b_base_j, lambda_j^2 * V_base_j)

b_target_j | eta_j, V_target_j
  ~ Normal(eta_j, V_target_j)

lambda_j ~ Uniform(0, 1)
```

The implementation uses the mean-field factorization
`q(eta_j, lambda_j) = q(eta_j) q(lambda_j)`. `q(eta_j)` is Gaussian and
`q(lambda_j)` is represented on bounded Gauss-Legendre quadrature nodes over
`(0,1)`. The optimized ELBO is

```text
E_q[log p(b_target | eta, V_target)
  + log p(eta | lambda, b_base, V_base)
  + log p(lambda)]
+ H[q(eta)] + H[q(lambda)].
```

Initialization sets the Gaussian mean and variance to the target Stage-1
posterior values; the first coordinate update derives `q(lambda)` from that
state. Updates stop when the absolute ELBO change is at most `m3.tol` (default
`1e-6`) or `m3.max_iter` (default `1000`) is reached. Log-domain normalization,
strictly interior quadrature nodes, and positive-variance validation provide
numerical safeguards.

An M3 chromosome is not checkpointed or scored if any calibrated variant has
not converged. Its diagnostic table is retained so the failed SNPs and iteration
counts can be inspected before retrying.

The returned calibrated effect is `BETA = E_q(eta_j)`. M3 is directional and
pairwise, does not read GWAS likelihoods, LD, annotations, phenotypes, or the
100 unselected Stage-1 candidates, and uses `VAR_BETA` directly without
squaring it.

### Final Stage 3

HERCULES scores the target validation and independent target test genotypes
with exactly:

```text
target_stage1_score    = X_target * b_target
calibrated_stage2_score = X_target * theta_target
```

The validation and test IIDs must be disjoint. The learner is fitted only on
target validation data, frozen, and then applied to target test predictors.
The optional test phenotype is used only for final R2/AUC evaluation.
Chromosome score files, the two predictor vectors, and phenotype tables must
cover the same cohort IIDs; HERCULES fails instead of silently shrinking a
cohort through an inner join.

- Quantitative: Lasso, ridge, and neural-network base learners.
- Binary: Lasso and neural-network base learners, binomial family, and
  `method.AUC` meta-learning.

## Platform and dependencies

Validated deployment target: Linux x86-64 with Python 3.11. Native Windows has
not been validated; Windows users should use WSL2, a Linux server, or a Linux
container/environment.

Required for the complete workflow:

- Python 3.11 and a C/C++ compiler for source installation;
- PLINK 1.9 and PLINK2;
- R and R packages `SuperLearner`, `glmnet`, `nnet`, `pROC`, and `data.table`;
- an OpenMP-capable compiler and BLAS are recommended.

## Installation

Clone or download the repository on Linux, then install from source:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install .
hercules --version
hercules doctor
```

Installing from source builds native extensions for the user's own Linux and
Python environment, avoiding dependence on a wheel built on another machine.
For development and tests:

```bash
python -m pip install -e '.[test]'
pytest
```

## Run the included examples

The repository includes deterministic synthetic data only; it contains no real
participants. From the repository root:

```bash
hercules config validate examples/data/hercules.quantitative.yaml
hercules run --config examples/data/hercules.quantitative.yaml
hercules run --config examples/data/hercules.binary.yaml
```

Edit the executable paths in the YAML files when PLINK, PLINK2, or Rscript are
not on `PATH`.

## Input files

### GWAS summary statistics

M1 and M2 require separate tab-delimited fastGWA-style files. Required columns:

| Column | Meaning |
|---|---|
| `CHR` | Chromosome |
| `SNP` | Unique variant identifier |
| `POS` | Base-pair position |
| `A1` | Effect allele |
| `A2` | Other allele |
| `N` | GWAS sample size |
| `AF1` | Effect-allele frequency |
| `BETA` | Marginal GWAS effect estimate |
| `SE` | Standard error |
| `P` | Ordinary GWAS association P value |
| `var_prior` | Optional positive, finite SNP-specific initialization variance |

Do not put a variance in `P`. HERCULES internally copies `var_prior` into the
initial variance field while retaining the ordinary P value. Missing
`var_prior` defaults to one. Duplicate SNPs, invalid alleles, and non-positive
or non-finite supplied variances are rejected.

Paths may contain `{chrom}` for chromosome-specific files.

### LD and genotype inputs

- `inputs.ld_reference.target`: target ancestry-matched magenpy LD store.
- `inputs.ld_reference.base`: base ancestry-matched magenpy LD store.
- `inputs.target_validation_genotype`: target validation PLINK prefix.
- `inputs.target_test_genotype`: independent target test PLINK prefix.
- `inputs.genotype_prefixes.base_validation` or
  `m2.validation_genotype`: base validation PLINK prefix.

Each PLINK prefix identifies `.bed`, `.bim`, and `.fam` files. The target
validation and test sample sets must not overlap.

### Phenotype inputs

`inputs.target_validation_phenotype` and optional
`inputs.target_test_phenotype` are tab-delimited tables with:

- `IID`;
- the column named by `inputs.phenotype_column`;
- every covariate listed in `inputs.covariates`.

The validation table trains Stage 3. The test table is never passed to model
fitting and is needed only to calculate the final metric.

M1/M2 candidate selection uses PLINK phenotype and keep files configured under
`m1` and `m2`. These are ancestry-matched validation samples, not the final
independent target test set.

### Annotation inputs

`inputs.functional_annotation` and `inputs.per_snp_heritability` record the
upstream annotation resources used to derive `var_prior`. The current public
runtime does not estimate SNP variances from raw annotation columns; the
positive derived values must be supplied in each ancestry's summary-statistics
file.

## YAML configuration

Start from [examples/hercules.example.yaml](examples/hercules.example.yaml).
Important fields are:

| Field | Description |
|---|---|
| `trait_name` | Safe output/run identifier |
| `chromosomes` | Chromosomes to process |
| `base_ancestry`, `target_ancestry` | Directional base and target labels |
| `inputs.summary_statistics.base_path` | M2 summary statistics |
| `inputs.summary_statistics.target_path` | M1 summary statistics |
| `inputs.ld_reference.base/target` | Ancestry-matched LD stores |
| `inputs.target_validation_genotype` | Target validation PLINK prefix |
| `inputs.target_validation_phenotype` | Stage-3 training phenotype table |
| `inputs.target_test_genotype` | Independent target test PLINK prefix |
| `inputs.target_test_phenotype` | Optional final-evaluation phenotype table |
| `inputs.phenotype_column` | Outcome column |
| `inputs.covariates` | Covariate columns |
| `inputs.trait_type` | `quantitative` or `binary` |
| `m1.validation_*` | Target Stage-1 model-selection files |
| `m2.validation_*` | Base Stage-1 model-selection files |
| `m3.max_iter`, `m3.tol` | Coordinate-optimization stopping settings |
| `m3.quadrature_points` | Bounded lambda quadrature resolution |
| `execution.seed` | Python/R deterministic seed |
| `tools.plink/plink2/rscript` | Executable names or absolute paths |
| `checkpoint.enabled/resume` | Versioned stage checkpoints |

The scientific settings `m3.model`, `m3.lambda_prior`, the 10 by 10 Stage-1
grid, and Stage-3 learner libraries are validated against the implemented
method and cannot silently select a different default model.

Configuration precedence is CLI override, environment variable, YAML value,
then package default where a command exposes an override.

## CLI

```bash
hercules --help
hercules --version
hercules doctor
hercules config validate hercules.yaml
hercules run --config hercules.yaml
hercules stage m1 --config hercules.yaml
hercules stage m2 --config hercules.yaml
hercules stage m3 --config hercules.yaml
hercules ensemble --config hercules.yaml
```

M3 requires completed selected M1/M2 posteriors. Ensemble requires M1 and M3
outputs. Checkpoints include the package version and scientific model IDs, so
checkpoints made by the replaced M3/ensemble implementation are not reused.

## Outputs

| Output | Contents |
|---|---|
| `HERCULES_M1.selected-posterior.tsv.gz` | Selected target Stage-1 posterior |
| `HERCULES_M2.selected-posterior.tsv.gz` | Selected base Stage-1 posterior |
| `HERCULES_M1/M2.selected-hyperparameters.tsv` | Selected grid settings |
| `HERCULES_M3.calibrated-posterior.tsv.gz` | One directional calibrated effect vector |
| `HERCULES_M3.convergence-diagnostics.tsv` | Lambda mean, convergence, iterations, ELBO |
| `HERCULES_ensemble.validation-scores.tsv` | Validation `IID`, p1, p2 |
| `HERCULES_ensemble.test-scores.tsv` | Test `IID`, p1, p2 |
| `HERCULES_ensemble.model.rds` | Frozen fitted Stage-3 model |
| `HERCULES_ensemble.predictions.tsv` | Independent test predictions |
| `HERCULES_ensemble.metrics.tsv` | Optional final R2 or AUC |

## R interface

Source [R/HERCULES.R](R/HERCULES.R) and call:

```r
HERCULES("/absolute/path/hercules.yaml")
```

This is a thin interface to the installed `hercules run` command. Python owns
configuration, stages, checkpoints, and external processes; R performs the
specified SuperLearner procedure.

## Validation status

The corrected workflow has automated unit tests and deterministic synthetic
Linux smoke tests. See [VALIDATION.md](VALIDATION.md) for commands and observed
results.

No authoritative historical result-generating Stage-2 outputs were available
for an old-versus-new manuscript-result comparison. Therefore:

> The implementation follows the scientific specification documented above,
> but equivalence to the historical manuscript results has not yet been
> demonstrated by numerical comparison.

## License and attribution

See [LICENSE](LICENSE), [NOTICE](NOTICE), and [CITATION.cff](CITATION.cff).
