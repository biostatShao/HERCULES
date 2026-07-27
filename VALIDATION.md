# Validation status

This file records checks that were actually executed for the clean source
distribution. It does not claim equivalence beyond the cases listed here.

## Validated environment

- Linux x86-64
- CPython 3.11
- source build through the PEP 517 build backend
- isolated installation outside the source checkout
- PLINK 1.9, PLINK2, and R 4.1.2 with the required ensemble packages
- deterministic example seed 7209, one inference thread, and one parallel job

Native Windows installation has not been validated. Windows users should use
WSL2, a Linux server, or a Linux container.

## Executed checks

- The source distribution and native wheel built successfully.
- The installed `hercules` package and all three native extensions imported.
- `hercules --version`, `hercules --help`, configuration validation, dependency
  diagnostics, individual stage commands, the ensemble command, and the full
  run command were exercised.
- All 57 automated tests passed.
- The quantitative and binary examples both completed end to end.
- The corrected quantitative workflow completed from a clean wheel
  installation outside the source checkout.

Observed example results:

| Example | Metric | Value | Wall time | Peak resident memory |
|---|---:|---:|---:|---:|
| Quantitative | R2 | 0.0577308741972629 | 68.94 s | 249,556 KiB |
| Binary | AUC | 0.683862433862434 | 19.32 s | 216,928 KiB |

These metrics are smoke-test outputs from synthetic data and are not estimates
of expected performance on real cohorts.

## Deterministic packaging comparison

The clean package was compared with the pre-unification reference workflow on
the same deterministic inputs. The comparison covered:

- M1 selected scores and all grid scores;
- M2 selected scores and all grid scores;
- M3 scores;
- M1, M2, and M3 posterior/fit tables;
- ensemble inputs, predictions, selected source, and metric table;
- row keys, sample identifiers, column names, shapes, and ordering.

All compared identifiers and schemas matched exactly. The maximum absolute
numeric difference was `0.0` for both the quantitative and binary workflows.
This demonstrates that the packaging and public-name unification did not alter
the deterministic example outputs.

That comparison predates the correction that makes the FastGWA `P` column a
fixed per-SNP prior variance throughout M1/M2 inference. It must not be cited as
an old-versus-new numerical equivalence result for the corrected prior-enabled
implementation.

## Corrected source-versus-wheel comparison

The corrected editable installation and a clean wheel installation were run on
the same quantitative fixture with seed 7209, one inference thread, and one
parallel job. The comparison covered:

- M1/M2 selected and grid scores;
- M1/M2 grid and selected posterior tables;
- selected hyperparameter and validation tables;
- M3 integrated posterior and scores;
- ensemble inputs, predictions, and metrics;
- all identifiers, columns, shapes, and ordering.

All compared schemas and non-numeric values matched exactly. The overall
maximum absolute numeric difference was `0.0`.

## Fixed per-SNP prior correction

The corrected workflow treats each stage's own FastGWA `P` column as a
precomputed per-SNP variance, converts it to prior precision using
`tau_beta = 1 / PVAL`, and prevents initialization or M-step updates from
replacing it. Input validation rejects non-numeric, non-finite, zero, and
negative values.

For the corrected quantitative and binary runs, the final mean `tau_beta`
reported by M1 and M2 was compared with `mean(1/P)` computed directly from each
stage's own FastGWA input. The absolute differences were `8.20e-7` for M1 and
`1.07e-6` for M2, consistent with float32 model/output precision. This confirms
that the configured prior reached the completed inference run and was not
replaced by the initialization or M-step update.

The recovered M3 implementation has no ancestry-bridging lambda variable and
no Beta/Uniform lambda update. No placeholder or newly invented lambda
algorithm has been added.

## Remaining scientific validation gap

One historical M3 comparison passed at `rtol=1e-8, atol=1e-10`. Historical M1
and M2 replays retained exact schemas and SNP order but did not pass the
predefined float32 tolerance. Raw annotation-to-per-SNP preprocessing also
remains an input-preparation responsibility and has not been reconstructed from
the original raw workflow.

Therefore:

> The refactor was designed to preserve the identified scientific defaults,
> but full scientific equivalence has not yet been demonstrated.
