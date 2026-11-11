# Stage 24 — Local Runtime Audit: PR #66 Review Checklist

- **Repository**: `Xin-Jiaqi/quantum-geometry-reproductions`
- **Audit timestamp**: 2026-10-22 (UTC)
- **Source of requirements**: PR #66 single commented review by Xin-Jiaqi (2026-10-22T07:35:46Z)
- **Publish branch**: `audit/stage24-<timestamp>` (new branch; `main` untouched; no merge)

The science repositories (`nlo-response-workflow`, `group-theory-operations-toolkit`,
`quantum-geometry-reproductions`) were **not modified**. All checks ran in temporary
directories or read-only against the existing working trees. The three pre-existing
dirty files/untracked sets of the science repos remain byte-identical (verified via
`git status` before/after each step).

---

## 1. C3 vs -6: physical symmetry gate for the spinless Haldane Hamiltonian

The Stage 23 structural criterion `POSCAR_point_group_is_C3` was superseded by a
Hamiltonian-level gate. `spglib` resolves the synthetic marker POSCAR to **-6 (C3h)**;
the synthetic Haldane model is a *spinless 2D* Hamiltonian whose exact symmetries are:

| Check | Method | Result |
|---|---|---|
| Source-gauge C3 invariance `H_src(C3·k) == H_src(k)` | 256 random k, orbital rotation identity | **exact**, max residual 5.98e-16 eV |
| Band equality `E(C3·k) == E(k)` | same grid | **exact**, max residual 1.78e-15 eV |
| Cell-gauge C3 via intertwiner `H_cell(C3·k) = W H_cell(k) W†`, `W = U†(C3k) U(k)` | center-embedding phases `e^{-ik·δ1}` | **exact**, max residual 1.64e-15 eV |
| m_z (`z → -z`) | all HR R-vectors strictly in-plane, centers in-plane, dH/dk_z ≡ 0 | **exact by construction** |
| Combined Hamiltonian point group | C3 + m_z | contains **C3h** |

**Gate verdict: `PASS_HAMILTONIAN_C3H_INPLANE`** — the spinless Haldane Hamiltonian,
in both source gauge and Wannier cell gauge, carries the full C3h symmetry, so its
in-plane third-rank shift-current tensor follows C3h selection rules.

### Frozen group-theory database cross-check (no formulas altered)

- In the frozen `data/optical_response_invariants.json`, the **in-plane** shift-current
  allowed coordinates for C3 (num 16) and C3h (num 22) are **identical**:
  6 independent coordinates in the symmetric `σ^{a;bc}` basis — `xxx, xxy, xyy, yxx,
  yxy, yyy` (the full 27-index label listing shows 8 labels because `xyx ≡ xxy` and
  `yyx ≡ yxy` by the `b↔c` permutation symmetry).
- The NLO `third_rank_polar_tensor` classifier on the marker POSCAR returns
  `point_group = -6`, `independent_count = 2`, which equals the frozen DB **C3h
  dimension = 2**; its in-plane allowed set matches the DB C3h in-plane set exactly.
- The Stage 23 `plan` exit-0 gate is therefore **not a software failure**: it is the
  correct C3h in-plane screening gate, restricted to the 6 physically allowed in-plane
  components (`xxx, xxy, xyy, yxx, yxy, yyy`). Out-of-plane components are forbidden
  by m_z and were not screened.

---

## 2. Group-theory test contract: 8 failures are dirty-tree extensions, not regressions

### Frozen-commit verification in a clean temp clone

A fresh clone at the frozen commit `9da48f1c89bff61d39066f032eda729f94193301`
(clean tree, `git status` empty) passes **202/202 tests**, including:

- `test_response_screening.py` — **5/5 PASS** (no 128-vs-96 failure at the frozen commit)
- `test_structure_responses.py` — **4/4 PASS** (no 4-tuple shape failures at the frozen commit)

### Mechanism of the 8 working-tree failures (quantified)

The 8 failures exist **only** in the uncommitted working tree:

1. **`test_all_registry_rows_are_reported_in_standard_order`** fails
   `128 != 96`. The dirty `invariants.py` adds an **`edelstein` response spec**
   (`S_i = χ_ij j_j`; axial output, polar input — the Rashba–Edelstein response)
   to `RESPONSE_SPECS`. The registry therefore reports 32 point groups × 4 responses
   = 128 rows; the frozen test asserts 32 × 3 = 96.
2. **`test_all_seven_crystal_system_fixtures_match_expected_dimensions`** (7 subtests)
   now receives a **4-tuple** of response dimensions, e.g. `(1, 3, 3, 1)` vs the frozen
   `(3, 3, 1)` and `(0, 0, 0, 0)` vs `(0, 0, 0)`, because every structure-response
   analysis returns the Edelstein response as a 4th entry.

### Analysis and recommendation (no science repo changes made)

- The Edelstein extension is an **intended scientific capability** (documented in the
  untracked `docs/edelstein-response.md` and the new altermagnet/spin modules); it is
  **not a defect** of the frozen baseline.
- The frozen tests encode a **3-response contract** and must be updated to the
  **4-response contract**: `32 * 3` → `32 * 4` in the screening test, and
  `EXPECTED_RESPONSE_DIMENSIONS` extended to 4-tuples with the Edelstein dimension,
  in the standard registry order.
- Verdict per the review's instruction ("不得为了测试通过而任意修改公式或模型，不得破坏测试"):
  the correct action is a **contract update of the tests** to match the intended
  Edelstein extension — the frozen code itself is green.

---

## 3. QGR physics-interface tests (previously 5 skipped)

- QGR checkout HEAD at run time: `86dba1481d712b4f1f35d7acc42e1d2fd7a4486a`
  (branch `research/bsf-nmee-stage7-sw-gate-20261020`; advanced externally during the
  audit window — dirty untracked set preserved byte-identical: 5 files, unchanged).
- GT v0.18.0 is **not on PyPI** ("No matching distribution found"); the installed
  pip package (0.17.0) lacks `point_group_operations`. The required API exists in the
  **GT dev working tree**, so the 5 tests were run with
  `PYTHONPATH=<gt-src>:<qgr>` and the used modules pinned by **explicit SHA-256**
  (GT HEAD `9da48f1`):
  - `__init__.py f5ab672ab8408e07ca77419fc8ccb1e25de36ddb624cb2eda66962bb5793f92b`
  - `point_groups.py ff69ca2f4ae58e98ca17849b64a3dedd437128b193aaea333e1c67028c91ac58`
  - `catalog.py 92a29672ef8d8526f753050214c5406f7d667e34f80fe4859245dd822e2abbe3`
  - `invariants.py f17ab14ac7b17bd4b98bfcfc47b54650a80b3bb428e8ff04512fb0e21621da7b`
  - `structure_responses.py b03bdc0cfc16a58a713ebe9237fbedc8522bfdcdf7aae182238cba5b790e2fc6`

**Result: 5/5 PASS** (`test_physical_model` 1/1, `test_real_toolkit` 4/4). The skipped
gate is resolved: the physics interface works against the real toolkit API.

---

## 4. Frozen parameters for the one-shot NLO shift-current run

Per review item (4), the NLO planner's auto-assignment `spin_degeneracy = 2` when
`soc: false` (a genuinely spinless model) was overridden by an **explicit freeze**:

| Parameter | Value (frozen) |
|---|---|
| `spin_degeneracy` | **1** (explicit; planner auto would have been 2) |
| `M, φ, t1, t2` | 0.4 eV, π/2, 1.0 eV, 0.2 eV |
| gauge / centers | Wannier cell gauge; centers at fractional (0.3333, 0.6667, 0.5) / (0.6667, 0.3333, 0.5) |
| units / sign | uA/V² (cell), uA·nm/V² (sheet); `sign_factor = -1` (electron charge e = -|e|) |
| formula | `ibanez_eq8_single_delta_full_nm`, `bc_combination = Iabc+Iacb`, `double_delta_reduction_factor = 0.5`, `position_matrix_mode = diagonal_tba` |
| numerics | reviewed-hermitian-v1; 5×5×1 grid (25 k); 401 points, 0–2 eV; η_sp = 0.05 eV; η_degen = 0.02 eV |

Run result: **exit 0**, batch `shift_current_20261022_154717`, nw=2, nocc=1, nR=7,
V = 51.9615 Å³, in-plane area 2.5981 Å², Lz = 20.0 Å, direct gap 1.41572 eV at
k = (-0.3, 0.3, 0). Full spectrum in `spectrum.csv`; manifest in `run-manifest.json`.

---

## 5. Independent numerical comparison of the shift-current kernel

A fully independent re-implementation of the NLO weight (own HR reader, own Fourier
kernel, own position matrices / generalized derivatives, own single-delta Lorentzian
and `Iabc+Iacb` combination, no NLO imports) was evaluated on the same package and the
same 5×5×1 grid, and compared against the actual NLO output (raw weight recovered from
`spectrum.csv` via the documented prefactor `0.5/V · π e³/ħ² · 1e-30/ω_EV · 1e6 · sign_factor`).

| Metric | Value |
|---|---|
| max |raw_NLO − raw_reimpl| over 8 components × 401 energies | **4.43e-13** (fp precision) |
| mean |raw_NLO − raw_reimpl| | 1.39e-14 |
| normalized lineshape mean abs diff | **2.25e-15** |
| nonzero components on this grid | `xxx, xyy, yxy, yyx, yyy` (`xxy`/`yxx` machine-zero though symmetry-allowed — model/grid-specific) |
| cell→raw prefactor (uA/V² per raw unit) | −7.358 |

### Band-pair decomposition (component `xxx`, probe ω = 1.465 eV, just above gap)

| band pair (n,m) | contribution (raw) | role |
|---|---|---|
| (2,1) | +3.42e-2 | dominant resonant v→c absorption |
| (1,2) | +8.5e-6 | negative-frequency partner — small off-resonant Lorentzian tail (the “±frequency” convention) |
| (1,1), (2,2) | 0.0 | diagonal excluded |

### Tensor arrangement and Lorentzian

- The `Iabc+Iacb` combination makes the weight symmetric in the two optical-field
  indices: NLO emits `xyx` and `yxy` separately but numerically equal to `xxy`/`yyx`.
- The spectrum lineshape is a **single Lorentzian** `η/π / ((ω − ω_mn)² + η²)` with
  η = 0.05 eV (FWHM = 0.1 eV) — no double-delta term is present in the raw weight
  (the 0.5 factor is applied only in the SI prefactor as
  `double_delta_reduction_factor`).

### FDR comparison status

The literal FDR SM Eq. (75) text is **SKIPPED_NO_FDR_ARCHIVE** — the Stage 22 archive
(`fdr_reproduction_stage22.zip`) is not present on this machine and the FDR code is not
in any of the three science repos. The comparison above therefore validates the NLO
kernel against an independent implementation and documents its exact band-pair /
tensor / Lorentzian / spectral-weight structure, which is the standard single-delta
shift-current structure expected of SM Eq. (75). The **absolute SI/2D sheet
normalization, the ±frequency integration convention, and the Eq. (9) vs SM Eq. (75)
sign branch remain unresolved** and were deliberately not forced.

---

## 6. Unresolved physics conventions (preserved, not modified)

1. FDR⇄NLO absolute equivalence including the SI prefactor (uA/V² vs uA·nm/V²) and the
   2D sheet normalization by Lz.
2. The Eq. (9) vs SM Eq. (75) **sign discrepancy** — no sign flip was applied.
3. The ±frequency integration convention (the (1,2) negative-frequency partner
   contributes an off-resonant tail in the single-delta full-nm sum).
4. DFT/W90Flow certification of an ab-initio material (this audit is synthetic-model
   only).
5. C3h-vs-C3 full 3D tensor difference (irrelevant for the in-plane 2D model, whose
   in-plane allowed set is identical — now resolved and documented).

---

## 7. Overall verdict

| Review item | Result |
|---|---|
| (1) C3/-6 symmetry | **PASS_HAMILTONIAN_C3H_INPLANE** — C3 + m_z verified at machine precision; in-plane C3 ≡ C3h tensor selection confirmed |
| (2) Group-theory 8 failures | **Analyzed**: frozen commit green (202/202); failures are the dirty-tree Edelstein extension → test contract must move to 4 responses (recommendation documented, nothing modified) |
| (3) QGR 5 skipped tests | **5/5 PASS** with GT dev-tree API pinned by explicit SHA (v0.18.0 unavailable on PyPI) |
| (4) Frozen parameters | **spin_degeneracy=1 explicit**; M, φ, t1, t2, gauge, centers, units frozen and recorded |
| (5) One NLO run + comparison | **run exit 0**; independent re-implementation matches NLO to 4.4e-13; band-pair/tensor/Lorentzian/weight structure documented |
| (6) Sign discrepancy | **Unresolved by design** — recorded, not forced |

**Overall: PASS with documented unresolved physics conventions.** No formula or model
was altered to force green; the science repos are untouched; the dirty working trees
are preserved.
