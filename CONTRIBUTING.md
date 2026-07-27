# Contributing

1. Do not change priors, grids, convergence, annotation weighting, M3 variance
   interpretation, ancestry integration, ensemble splitting, or binary
   selection without explicit scientific authorization.
2. Add a differential test before changing numerical code.
3. Engineering changes should preserve schemas, ordering, and defaults.
4. Numerical changes require before/after runtime, memory, and explicit
   tolerances.
5. Keep subprocesses as argument vectors, check return codes, and retain logs.
6. Do not add private data, credentials, machine-specific defaults, generated
   binaries, or local build outputs.
7. Run `pytest`, `python -m build`, clean-install imports, naming/path searches,
   and the relevant real workflow smoke test before submitting a change.
8. Preserve all copyright and attribution notices in `LICENSE` and `NOTICE`.

