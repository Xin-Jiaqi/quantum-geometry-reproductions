# Stage 23 — Local Runtime Audit

**Overall status:** `FAIL_REVIEW_REQUIRED`

- Python: `/Library/Frameworks/Python.framework/Versions/3.12/bin/python3` (3.12.4)
- Report JSON: `<workspace>/.stage23_work/stage23_20261022-142840/results/STAGE23_LOCAL_RUNTIME_AUDIT.json`
- Logs: `<workspace>/.stage23_work/stage23_20261022-142840/results/logs`

## Command results

| Check | Status | Duration (s) | Return code |
|---|---|---:|---:|
| nlo_unittest_discover | `PASS` | 1.052 | 0 |
| nlo_actual_cli_preflight | `FAIL_OR_BLOCKED` |  |  |
| group_theory_unittest_discover | `FAIL` | 136.288 | 1 |
| group_theory_database_validate | `PASS` | 0.606 | 0 |
| qgr_physics_interface_registry_tests | `FAIL` | 0.062 | 1 |
| fdr_stage22_read_only_pilot | `SKIPPED_NO_FDR_ZIP` |  |  |
| fdr_full_verify_py | `SKIPPED_NO_FDR_ZIP` |  |  |

## Repository integrity

| Repository | Commit before | Commit after | Worktree unchanged |
|---|---|---|---|
| nlo | `dc9eec064c68e103eb974af7bc79eea482c05dfb` | `dc9eec064c68e103eb974af7bc79eea482c05dfb` | `True` |
| group_theory | `9da48f1c89bff61d39066f032eda729f94193301` | `9da48f1c89bff61d39066f032eda729f94193301` | `True` |
| qgr | `682a1591854a1bc38f4102ae5b0f056fc84b364c` | `682a1591854a1bc38f4102ae5b0f056fc84b364c` | `True` |

## Scientific boundary

The NLO integration fixture is a synthetic two-band Haldane tight-binding model with an exact HR Fourier mapping; it is not an ab-initio W90Flow export. The audit does not claim equivalence of the FDR and NLO response formulas or absolute SI normalization. No NLO response `run` command is executed.

---

## Supplementary evidence (appended by the Stage 23 runner)

### Test counts (from real command output)

| Suite | Total | Passed | Failed | Skipped | Errors | Duration (s) | Exit |
|---|---:|---:|---:|---:|---:|---:|---:|
| NLO `unittest discover -s tests` | 31 | 31 | 0 | 0 | 0 | 1.052 | 0 |
| Group-theory `unittest discover -s tests` | 243 | 235 | 8 | 0 | 0 | 136.288 | 1 |
| Group-theory database validate | 1 | 1 | 0 | 0 | 0 | 0.606 | 0 |
| QGR physics-interface registry (temp clone, branch `integration/physics-capabilities-handoff-v1`) | 12 | 12 | 0 | 0 | 0 | 0.004 | 0 |
| QGR `test_physical_model.py` (current checkout) | 1 | 0 | 0 | 1 | 0 | <1 | 0 |
| QGR `test_real_toolkit.py` (current checkout) | 4 | 0 | 0 | 4 | 0 | <1 | 0 |
| FDR Stage 22 pilot / `verify.py` | — | — | — | — | — | — | SKIPPED_NO_FDR_ZIP |

### Haldane → HR mapping and NLO kernel (real numbers)

- Standalone exact synthetic mapping: **11/11 checks pass** (128 sampled k points; max residuals: H 9.04e-16 eV, eigenvalues 1.78e-15 eV, dH 9.99e-16 eV·Å, d2H 8.88e-16 eV·Å²; HR round-trip and real-space Hermiticity exactly 0; min sampled direct gap 1.3723 eV).
- NLO HR reader/Fourier kernel on the same package (64 sampled k points): **7/8 checks pass** — H reconstruction 1.02e-15 eV, eigenvalues 1.33e-15 eV, dH 8.77e-13 eV·Å, d2H 1.62e-12 eV·Å², raw HR pairing 0, gap open 1.2908 eV.
- The 8th NLO check **fails**: `POSCAR_point_group_is_C3 = false`. spglib 2.5.0 resolves the synthetic marker POSCAR to point group **-6** (6 operations), not the intended C3 (`3`). Per the kit's own rule, the NLO `plan` must be treated as **blocked** until this symmetry fixture is corrected. This is a real, reproducible finding of this environment.

### NLO `inspect` and `plan` (real CLI output, corrected invocation)

| Command | Exit | Key output |
|---|---|---|
| `nlo-response inspect <package> --structure L1` | 0 | num_wann=2, nocc_model=1, dimensionality=2D, `input_integrity.hashes_verified=true` (SHA256 match against package-manifest.json), `symmetry.point_group="-6"`, 8 allowed nonzero components (`xxx`,`xyx` independent) |
| `nlo-response plan <package> --structure L1 --response shift_current --kmesh 3x3x1 ... --no-mat` | 0 | `config.response_type=shift_current`, `nk=[3,3,1]`, 8 tensor components, `position_matrix_mode=diagonal_tba`, no spectrum computed |

> The first audit run failed these two commands because the kit script passed a spurious literal argument `nlo-response` to `python -c` (argparse: `invalid choice: 'nlo-response'`). Per the execution prompt's rule 3 the kit script itself was fixed (removed the token) and the audit was re-run; no science repository was touched.

### QGR physics-interface tests — honest status

- `integration/physics_interface/tests` does **not exist** in the user's checked-out QGR branch (`research/bsf-nmee-stage7-sw-gate-20261020`); the kit-targeted path exists only on remote branch `integration/physics-capabilities-handoff-v1` (tip `73e0686`). The kit audit therefore recorded `FAIL` (exit 1, `ImportError: Start directory is not importable`).
- Supplementary: in a **separate temporary clone** at that integration branch, the registry suite ran **12/12 PASS** (protocol registry integrity, no physics-dispatch claims, safe live files/evidence, source-revision requirements).
- The physics-interface tests present in the current checkout (`test_physical_model.py`, `test_real_toolkit.py`) are **skipped** (5/5) because the installed `group-theory-operations-toolkit` is 0.17.0 while the tests require 0.18.0 / `point_group_operations`. No dependency was installed.

### Group-theory failures (not caused by Stage 23)

8 real failures in the repository's working tree (which already carried 16 pre-existing uncommitted changes before this audit):
- `test_response_screening...test_all_registry_rows_are_reported_in_standard_order`
- `test_structure_responses...test_all_seven_crystal_system_fixtures_match_expected_dimensions` on 7 of 7 crystal-system fixtures — `AssertionError` shape-convention mismatch, e.g. `(1,3,3,1)` vs expected `(3,3,1)` (4-tuple vs 3-tuple response tensor dimensions).

The group-theory database validation (`load_database` + `validate_database`) reports **no errors** (exit 0).

### Repository integrity after the whole audit

All three repositories keep their original HEAD, branch, and dirty state (NLO 0→0, group-theory 16→16, QGR 3→3). No checkout, pull, fetch, reset, stash, commit, or push was performed on any of the three science repositories.

### Stage 22 FDR acceptance

`fdr_reproduction_stage22.zip` was **not found** on this machine (searched Downloads, Desktop, Documents, and the home tree). The FDR read-only pilot and `verify.py` regression are recorded as `SKIPPED_NO_FDR_ZIP`. They must be run after the archive is provided.

### Scientific judgement (unchanged scope)

This run verifies the mapping *解析 Haldane Hamiltonian → HR file → NLO Fourier kernel* for a synthetic two-band model (H, bands and position-dressed first/second k-derivatives all match to ≤1.6e-12 in the stated units). It does **not** establish that the FDR shift-current formula equals the NLO diagonal-TBA shift-current kernel: the absolute SI prefactor, charge/sign convention, 2D sheet normalization, frequency integration factors, broadening protocol, tensor ordering, and the main-text Eq. (9) vs SM Eq. (75) sign discrepancy all remain unverified. The synthetic fixture's intended C3 symmetry is not realized under the local spglib, so the NLO plan gate is blocked until that fixture is corrected.
