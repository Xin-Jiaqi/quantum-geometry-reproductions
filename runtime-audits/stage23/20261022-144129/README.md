# Stage 23 Local Runtime Audit — 20261022-144129

This directory contains the machine-generated evidence for the Stage 23 local runtime
acceptance of the three research repositories, produced on a local macOS (arm64)
environment. All files here were sanitized for publishing (home paths, temp paths and
host names replaced; no credentials, tokens or private keys are present).

## Test environment

- OS: Darwin 27.0.0 (arm64), Python 3.12.4 (framework build)
- Packages: numpy 1.26.4, scipy 1.14.0, spglib 2.5.0, ase 3.23.0,
  materials-structure-core 0.0.2, pytest 9.1.1;
  pip-installed `group-theory-operations-toolkit` 0.17.0 (the QGR physics tests demand 0.18.0)
- No package was installed, no network was required, no science repository was modified.

## Repositories audited (local clones; HEADs at audit time)

| Repository | Local checkout | Branch | HEAD |
|---|---|---|---|
| nlo-response-workflow | <git-checkouts>/nlo-response-workflow | codex/research-sync-20261019 | dc9eec064c68e103eb974af7bc79eea482c05dfb |
| group-theory-operations-toolkit | <git-checkouts>/group-theory-operations-toolkit | main | 9da48f1c89bff61d39066f032eda729f94193301 |
| quantum-geometry-reproductions | <workspace>/quantum-geometry-reproductions | research/bsf-nmee-stage7-sw-gate-20261020 | 682a1591854a1bc38f4102ae5b0f056fc84b364c |

All three were unchanged (HEAD, branch, dirty state) before and after the audit.

## Acceptance scope and results (summary)

| Item | Result |
|---|---|
| Stage 23 kit bridge tests (3) | 3/3 passed |
| Haldane → Wannier90 HR exact mapping (128 k points) | 11/11 checks passed (synthetic model) |
| NLO HR reader + Fourier kernel (64 k points) | 7/8 checks passed; point-group C3 check FAILED (observed `-6`) |
| NLO `inspect` | exit 0, hashes verified, point group `-6`, 2D |
| NLO `plan` (shift_current, 3x3x1, no-mat) | exit 0 (pipeline plumbing only; no `run`) |
| NLO repository `unittest` | 31/31 passed |
| group-theory repository `unittest` | 235/243 passed, 8 failed (shape-convention assertions) |
| group-theory database validation | passed (no errors) |
| QGR physics-interface registry tests | 12/12 passed in a separate temp clone at branch `integration/physics-capabilities-handoff-v1` (the path is absent from the audited checkout branch); current-checkout physics tests skipped (require toolkit v0.18.0) |
| Stage 22 FDR pilot + `verify.py` | SKIPPED — `fdr_reproduction_stage22.zip` not present on this machine |

The NLO plan gate is intentionally reported as **blocked** because the synthetic
marker POSCAR resolves to point group `-6` under the local spglib, not the intended C3.
Per the kit's own rule the fixture must be corrected before plan-based runs are meaningful.

## Files

- `STAGE23_LOCAL_RUNTIME_AUDIT.md` / `.json` — kit-generated audit report (JSON contains
  per-command exit codes, durations, repo before/after snapshots; MD is the summary).
- `STAGE23_SUPPLEMENTARY_EVIDENCE.json` — runner-added counts, QGR supplementary runs,
  environment, and the scientific boundary.
- `STAGE23_HALDANE_WANNIER_MAPPING.json` — the exact synthetic Haldane→HR mapping record.
- `logs/` — sanitized per-command stdout/stderr logs.
- `models/stage23_haldane_synthetic/` — the small synthetic two-band Haldane fixture
  (HR, POSCAR, w-center, quality, manifest). It is a tight-binding test fixture only,
  NOT an ab-initio/W90Flow material model.
- `SHA256SUMS.txt` — SHA-256 checksums of every file in this directory.

## Limitations (what this audit does NOT claim)

- Interface/pipeline passes do not establish equality of the FDR and NLO response formulas.
- Absolute SI prefactor, charge/sign convention, 2D sheet normalization, frequency
  integration, broadening protocol and tensor ordering are unverified.
- The main-text Eq. (9) vs SM Eq. (75) sign discrepancy remains open.
- No DFT/W90Flow Wannierization is certified; no NLO response `run` was executed.
- Stage 22 FDR regression could not be run (archive missing on this machine).
