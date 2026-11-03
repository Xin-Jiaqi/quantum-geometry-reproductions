# MATERIAL_HSA_SYNTHETIC_IMPORT_AUDIT — 2026-10-23

Verdict: **PASS_SYNTHETIC_FORMAT_PARSER** (NOT `MATERIAL_AUTHENTICATED`)

Machine-readable receipt: `MATERIAL_HSA_SYNTHETIC_IMPORT_AUDIT.json` (schema
`hsa-synthetic-import-audit-v1`, written by `scripts/test_import_hsa_operator_bundle.py`).
All work used only ordinary Python/NumPy on the existing four-band synthetic
model; no VASP/PBS/Wannier executable or new AgF₂ computation was run
(`no_wannier_executed: true`). Local CPU ≈ 4.5 min for the import audit.

---

## 1. Phase A — synthetic JSON carrier roundtrip (SYNTHETIC_CARRIER_PASS)

Carriers written with a manifest `file_sha256` (canonical `sort_keys`,
`separators=(',',':')` bytes hash) for `build_v2()` and
`build_v2(position_offdiag=0.3·np.kron(SX,SZ))` carrying
H(R)/S(R)/A(R)/lattice/R/centers/units eV/Å with `spin_axis=(0,0,1)` and
`position_center_once=true`. The strict reader verifies: hash, shapes,
finite values, R↔−R pairing Hermiticity, spin Hermiticity and spectral norm
≤ 1+1e-9, and **centers == diag A(R=0) counted exactly once** (R=(0,0) need
not be the first R entry).

- Roundtrip recovery: **0.0** for all three operators on both fixtures
  (bitwise through canonical JSON).
- Responses at 41² from the re-read carriers:
  - v2_positive direct **+0.0342937475**, Wilson slope **+0.0341517445**
    (both within 1e-10 of the reference anchors).
  - v2_spinZ_pos direct −0.02536 / Wilson −0.02557 (self-consistent with its
    source fixture; no reference is claimed for this NC-style fixture).
- Negative controls — every one must raise (nonzero exit), and all did:

| negative | detected by | observed |
|---|---|---|
| corrupt hash | manifest sha256 mismatch (checked before parsing) | hash mismatch |
| double centers | centers vs diag A(0) | ValueError |
| remove position_r | required matrix missing | ValueError |
| wrong units (meV) | units dict equality | ValueError |
| non-Hermitian spin | R/−R Hermiticity of S | ValueError |
| missing −R partner | R-list structural check before matrix parse | ValueError |
| spin-gauge swap (S rotated alone) | response | 0.03429→0.02141 (Δ=0.0129, 38 %) |
| wrong lattice orientation (a1/a2 rows swapped) | response | →1.2e-17 (Δ=0.0343) |

## 2. Phase B — genuine Wannier90 syntax, synthetic files (SYNTHETIC_FORMAT_PARSER_PASS)

`scripts/wannier90_synthetic_fixtures.py` authors official-syntax
hr.dat/r.dat/spn/u.mat/u_dis.mat/v.mat/eig/win files for the supported case
(2D slab, orthorhombic, Rz=0, `write_ndegen_applied = .true.` deg=1, no WS
fold, `saxis = 0 0 1`, spinors, single explicit k ordering, formatted spn
with the documented `m,n,i` loop order, i fastest). The parser
`reconstruct_bundle_from_w90` builds V from v.mat (V=U_dis·U) > u_dis.mat >
u.mat, checks V†V−I ≤ 1e-8 and H(k) vs V†diag(ε)V at every k, forms
S^W(k)=V†S^B(k)V(k) and r.dat **directly as the final gauge** (no second
U transform; centers = diag A(R=0)).

- **B1 identity-physics gauge at 41²** (final operators == the source
  bundle): H/A recovery 0.0, S 1.0e-10 (16-decimal spn-file floor),
  H-vs-eig 2.3e-12, direct 0.0342937475, Wilson 0.0341517445.
- **B2 smooth nonidentity gauge at 21²** (361-R window; k-dependent W with
  the exact analytic position-connection term A'=W†AW + iW†∂W): H-vs-eig
  7.7e-10, V†V−I 4.8e-15, direct 0.0342937475 at both 21² and 41² — the
  joint k-dependent gauge incl. the connection law **preserves the
  response**; Wilson 21² 0.0337756 vs source 0.0337753 (documented
  mesh-truncated drift, |Δ|=2.1e-5).
- **B3 disentangled 6-band at 21²** (V=[U_eig†W;0], U_dis=[U_eig†;0],
  ε + far bands 25/26): H-vs-eig 5.3e-10, V†V−I 3.3e-15, direct
  0.0342937475.
- Negatives (fail closed):
  - wrong U orientation (per-k transposed u.mat): raises
    "H(k) mismatch vs eig" at k=0 (err 0.66 ≫ tol).
  - mismatched k order (spn k-blocks permuted): reconstructs, but the
    response shifts by 8.1e-6 (5 orders above the 1e-11 numerical floor) —
    a clear, documented weak signature of the wrong k pairing.

Tolerances are recorded in the JSON receipt: B1 op recovery ≤ 2e-10,
B2/B3 H-vs-eig ≤ 1e-8 at the 361-R window, B2 direct ≤ 1e-5 vs the 41²
anchor, B2 Wilson ≤ 1e-3 vs the same-mesh source.

## 3. Phase C — Wilson winding fix (WILSON_WINDING_FIXED)

The old predicate `abs(wrap(last−first)) >= pi` is unreachable because
`wrap ∈ [−π, π)`. Replaced by `wilson_winding_unwrapped`: integer
`round((u[-1]−u[0])/2π)` of `np.unwrap(rows)` over the periodic ky cycle,
with closure residual and a jump-like flag; an in-range boundary jump
(principal branch of γ(0)+2πC) is resolved by re-reading the interior
drift.

- Regression 1 (exact): `γ = 2π·ky` incl. the closing link — winding +1,
  closure residual 0.0, not jump-like; the old predicate cannot fire.
- Regression 2 (physical): 2-band Chern model (m=0.5, C=−1) Wilson rows —
  winding ±1 recovered via the interior drift, closure residual 0.27 rad
  (mesh-truncated), not jump-like.
- v2_positive: winding 0, closure 0.0017 rad, gate's `wilson_branch_check`
  passes (no spurious Chern flag).

`test_multiband_igmc_gate.py`'s `wilson_branch_check` now records per-b
`winding_integer / winding_closure_residual_rad / winding_max_local_step_rad
/ winding_jump_like` and sets `chern_winding_flagged` from any nonzero
winding; `test_gate_regressions.py` R4 comment corrected ("RAW drift 3.5 rad
> pi (wrapped ~2.78 rad < pi: the raw/unwrapped difference is the
branch-jump-in-b detector); rows off cut").

## 4. Honest limitations

- `MATERIAL_AUTHENTICATED` is **never** claimed: the parser is validated on
  synthetic authoring only; upstream has LWAVE=.FALSE. so no WAVECAR/W90
  files exist to import, and none are asserted.
- Smooth/disentangled final gauges need the 361-R window because the
  k-dependent rotation's Fourier content is not confined to the source's
  5-R list; the response invariance is exact for the identity gauge and
  ≤ 1e-10 (direct) / 2.1e-5 (Wilson) for the smooth gauge at its window.
- The spn-k-permutation negative is weak by design of the smooth spin model
  (8.1e-6 ≫ 1e-11 floor); a material spin texture would show a much larger
  shift, but no material claim is made here.

## 5. Files

- `scripts/import_hsa_operator_bundle.py` — carrier + W90 parsers +
  `wilson_winding_unwrapped` (Phase C).
- `scripts/wannier90_synthetic_fixtures.py` — official-syntax authoring.
- `scripts/test_import_hsa_operator_bundle.py` — phases A/B/C audit.
- `scripts/test_multiband_igmc_gate.py`, `scripts/test_gate_regressions.py`
  — gate regression (re-run after the branch-check change).
- `scripts/wannier_operator_contract_pinned26.py` — unchanged pinned module
  (sha256 asserted by R6 in the gate regressions).

## 6. Gate regression after the branch-check change (PASS)

`python3 scripts/test_gate_regressions.py` → **7/7 passed** (R1 wrong-factor
NC4 guard, R3 branch-unstable/cut margin, R4 finite-b raw-drift detector
with the corrected comment, R5 dirac_2band anchors 0.0152366470 /
0.0151827707, R6 pinned-module sha256, R7 direct21/Wilson21 + corrected NC4
zero).
`python3 scripts/test_multiband_igmc_gate.py` → **PASS** (236 s): direct
+0.0342937475 vs Wilson +0.0341517445 at 41² (0.41%), gauge invariance
residual 2.8e-17 / 1.1e-12, NC1–NC5 all detected, `wilson_branch` now
records per-b `winding_integer` (all 0), closure residuals ~1.7e-3 rad,
max local steps ~0.022 rad, `chern_winding_flagged: false` — the v2_positive
Wilson rows are smooth and the winding fix does not disturb the hardened
gate.

## 7. Next step (unique recommendation)

The import path is now validated end-to-end on synthetic official-syntax
files with the four-band model; the AgF₂ H/S/A import itself remains blocked
on an authorized real-material Wannier90 run (upstream has LWAVE=.FALSE., so
no WAVECAR/W90 files exist to import; out of scope of this round). The single
next step is to close PR #38 with this audit as the import-path deliverable
(verdict PASS_SYNTHETIC_FORMAT_PARSER, material_authenticated=false) and
defer any real-material import until an authorized DFT/Wannier round
provides authentic hr.dat/r.dat/spn/u.mat/v.mat from a Wannier90 execution.
