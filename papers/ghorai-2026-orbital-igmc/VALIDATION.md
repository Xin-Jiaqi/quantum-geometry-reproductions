# Local validation record

2026-11-09; Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0.

`PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v`:
**23 tests passed, no skips** (21 internal + 2 cross-paper bridges).
`reproduce_core.py` was executed; numerical arrays are stored in `results/core.json`.

Cross-paper files were fetched from baseline commit `9ffa75d3d48d8ff9c2152a43eebfb3cbf705ff5a` and verified byte-for-byte using Git blob hashes before the bridge tests:

- NMEE `nmee_lab/ptbackend/cumnas.py`: `145d756e01ed4f63268683f8960c5031255cd143`.
- Spin IGMC `src/igmc_lab/core.py`: `804bb77c66b23049b296f4eda64782cfdcb37916`.

The local checkout contains the new module and the two reference files, not the entire repository. This is **not** a claim that all six older suites were rerun. Python 3.11/3.12 full-repository jobs are configured in the GitHub workflow; their actual status must be read from the PR, not inferred from this local run.

The tests certify explicit numerical protocols and the declared symmetry/invariance checks. The orbital D10/A9 discrepancy remains open; passing tests do not settle which absolute physical amplitude is correct.
