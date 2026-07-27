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
- All 48 automated tests passed.
- The quantitative and binary examples both completed end to end.

Observed example results:

| Example | Metric | Value | Wall time | Peak resident memory |
|---|---:|---:|---:|---:|
| Quantitative | R2 | 0.019599803263548 | 35.76 s | 247,732 KiB |
| Binary | AUC | 0.683862433862434 | 26.39 s | 216,956 KiB |

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

## Remaining scientific validation gap

One historical M3 comparison passed at `rtol=1e-8, atol=1e-10`. Historical M1
and M2 replays retained exact schemas and SNP order but did not pass the
predefined float32 tolerance. Raw annotation-to-per-SNP preprocessing also
remains an input-preparation responsibility and has not been reconstructed from
the original raw workflow.

Therefore:

> The refactor was designed to preserve the identified scientific defaults,
> but full scientific equivalence has not yet been demonstrated.
