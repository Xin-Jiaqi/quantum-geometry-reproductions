# Round 9 worklog | 2026-10-20

## Scope

Publication-preparation of independently validated clean two-band IGMC calculations. No additional empirical parameter fit or modification to preceding source, model, theory derivation or frozen numerical values.

## Local checkpoint

- Root: `/mnt/data/igmc_reproduction_round9`.
- Paper module: `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current`.
- The prior `scripts/verify_release.py` still calls all Round1–8 physics gates and Python tests. The strict Round8 frozen source/evidence SHA gate passed for 26 existing file paths (new R9 files are additive only).
- Added `release/release_contract.json`, `scripts/verify_round9.py`, `tests/test_round9_release.py`, `docs/RELEASE_CANDIDATE_ROUND9.md`, and staged `.github/workflows/igmc.yml`.
- `python scripts/verify_round9.py`: `PASS_LOCAL_DRAFT`, not GitHub CI. `python -m unittest discover -s tests -v`: **76 discovered, 75 passed, 1 skipped** (native toolkit optional).
- No APS article or SM PDF, original source scans, original user notes, or private source archives placed in the proposed public paper-module directory.

## Preserved physics

- Eq.12 in-gap Dirac D: `0.01775389 A/(m T)` at cyclic `f=1 THz, g=2`, not author Fig.2(d) raw amplitude.
- Eq.14 time-reversal response convention remains `PARTIAL`; printed Table II `-42'm'` equality remains source-table open.
- SM intermediate Eq.34/46 statements remain separate provenance questions, never inserted as empirical plot corrections.

## Next gate

1. Review the candidate public file tree and repository-code license choice.
2. Only with an explicit user go-ahead submit the new directory and new workflow on a branch/PR, leaving frozen NMEE untouched.
3. Record GitHub commit SHA and actual Actions result. Without step 2–3, publication and CI must remain `NOT DONE`.
4. Optionally run the pinned full group-theory toolkit API in a network-enabled environment to remove the one SKIP.

## Packaging qualification

- Source-only, public-staging-shaped ZIP includes 126 files and `.github/workflows/igmc.yml`; no publisher PDF/SM scan or private external `USER_PROTOCOL_SNAPSHOT`.
- ZIP extracted under `/mnt/data/_igmc_round9_zip_verify/`; `verify_round9.py` passed.
- Fresh Python 3.13 venv initially lacked setuptools (`BackendUnavailable: Cannot import setuptools.build_meta`). Offline resolution provided host-preinstalled build/dependency site packages explicitly via `PYTHONPATH=/opt/pyvenv/lib/python3.13/site-packages`, then installed the extracted package editable with `--no-build-isolation --no-deps`. Both import and release verification succeeded. This does **not** count as dependency-complete network-free cleanroom installation.
