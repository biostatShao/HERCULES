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

## Fixed per-SNP variance prior for M1 and M2

Both M1 and M2 use a fixed, stage-specific prior variance for every SNP. For
FastGWA input, the file's `P` column is deliberately repurposed: it must contain
the precomputed per-SNP variance (for example, the `baseline` value from an
LDSC-derived `.snpvg` file), not a statistical P-value. `magenpy` exposes this
input column internally as `PVAL`, and HERCULES computes

```text
tau_beta_j = 1 / PVAL_j
```

because `tau_beta` is the prior precision. The precision is fixed through
initialization and every M-step. M1 reads the target-ancestry FastGWA file; M2
independently reads the base-ancestry FastGWA file. Values must be finite and
strictly greater than zero.

```yaml
m1:
  per_snp_prior:
    enabled: true
    source: summary_statistics
    column: PVAL
    input_type: variance
    fixed_during_inference: true
m2:
  per_snp_prior:
    enabled: true
    source: summary_statistics
    column: PVAL
    input_type: variance
    fixed_during_inference: true
```

Do not place ordinary association P-values in this column when running
HERCULES.

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
The corrected package passed 57 automated tests. Quantitative and binary
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

The final reported mean `tau_beta` matched `mean(1/P)` from the stage-specific
input to float32 output precision: absolute differences were approximately
`8.2e-7` for M1 and `1.1e-6` for M2.

One historical M3 comparison passed at `rtol=1e-8, atol=1e-10`. Historical M1
and M2 replays retained exact schemas and SNP order but did not pass the
predefined float32 tolerance. Raw annotation-to-per-SNP preprocessing remains
an upstream input-preparation responsibility. The example FastGWA files use
positive synthetic per-SNP variances in their `P` columns.

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
