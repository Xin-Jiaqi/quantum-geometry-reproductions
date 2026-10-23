# Multiband same-gauge IGMC gate: hardening round (zero VASP)

**2026-10-23, task `HARNESS_SYNTHETIC_GATE_HARDENING_20261023.md`. Verdict: PASS (conditional semantics).** The corrected multiband low-omega Dxx result **survives** (direct +0.0342937 vs A-embedded Wilson slope +0.0341517 at 41², relative 0.41%), with all four hardening requirements met. PASS means the gate may now be used to develop a **genuine source-authenticated material H/S/A file importer**; it is **not** material IGMC, and **not** authorization for costly DFT/Wannier. No VASP/PBS/Wannier/SCF was executed. CPU: full gate ~176-193 s, regressions ~10 s (≤ 15 min total budget).

Repository: PR #38, branch `research/agf2-b1-wsl-relax-20261021`, deliverables at `idea_screen_1/soc_static_1/`.

## 1. What was fixed (each with evidence)

### 1.1 NC4: the old 'no-SOC, spin-conserving' control was not spin-conserving (mandatory)

Old `build_no_soc()` used `H = I⊗h2(k)`, `S = I⊗σ`: h2 and σ acted on the **same tensor factor**, so `||[H(k), S^x]||₂ = 2β = 1.2 eV` at k=0 (independently reproduced; the researcher's `scripts/audit_source_gate_invariants.py` runs clean against this branch and prints exactly `source_nc4_commutator_norm_eV = 1.2`). The old D≈0 was therefore an accidental cancellation, not evidence of a conserved-spin zero.

Fix: the gate now uses `build_no_soc_corrected()` with `H = h2(k)⊗I_spin`, `S_j = I_orb⊗σ_j`. Commutators asserted < 1e-10 eV for **all three** S_j at a 21² grid (computed max = 0.0). Re-computed on the corrected fixture at 41²:

- direct Dxx/ω = **1.1e-33** (machine zero), finite-field Wilson slope = **0.0** — both vanish, now **with** the conservation proof (cross-spin velocity elements and same-spin σ^x elements vanish in the occ–empty pair sum).
- The **wrong_factor** builder is retained as `build_no_soc()` solely for the injected-fault regression (R1): the gate's conservation assertion must reject it.

### 1.2 NC2: value-to-label inversion (mandatory)

The committed JSON/source are authoritative: for **v2_spinZ_pos**, `with_A (covariant) = −0.025357`, `omit_A (bare) = +0.034294`, |Δ| = 0.059651. The commit-adcd75e Markdown table reversed these labels. Repaired in this round's `MULTIBAND_SAME_GAUGE_IGMC_GATE.md` (dated erratum) and re-keyed in the JSON as explicit `with_A_covariant` / `omit_A_bare` while retaining the old keys. The physically correct operator set includes A; its response on the spin-Z position fixture is −0.02536. `v2_positive` (a **different** position connection) keeps +0.03429 with and without A (documented model insensitivity of the Â channel there); the response-level omit-A detection therefore lives on `v2_spinZ_pos`.

### 1.3 Global indirect gap (mandatory)

`gap_scan` now computes, per field b ∈ {0, ±1e-4, ±1e-3} and per mesh 21²/31²/41², the **true global indirect gap** `min_k ε₃(k) − max_k ε₂(k)` and **raises** (`global insulator FAILED`) when it is ≤ 0. Values (41²):

| fixture | indirect gap at b=0 | at ±1e-3 |
| --- | --- | --- |
| pinned_round26 | 1.0225 eV | 1.0225 eV |
| v2_positive | 0.9853 eV | 0.9853 eV |
| v2_spinZ_pos | 0.9853 eV | 0.9853 eV |
| no_soc_control (corrected) | 1.2000 eV | 1.1980 eV |

(21² gives 0.9892 eV for v2_positive; mesh-converging to 0.9853 at 41² — consistent with the helper's 21² diagnostics.) The injected-fault regression R2 builds a fixture with **positive direct gap at every k but negative indirect gap** (−0.28 eV) and requires rejection.

### 1.4 Embedded Wilson phase branch (mandatory)

`wilson_branch_check` now examines the **embedded** row phases actually used for Dxx/ω (A-transported rank-2 det-W rows) at b ∈ {0, ±1e-4, ±1e-3}:

- adjacent-ky principal-branch jumps: max 0.0225 rad at 41² (< π);
- branch margin: all row phases ≥ 0.94 rad from the ±π cut (branch-unstable rows within 0.05 rad raise);
- finite-b continuity (raw, unwrapped drift of γ(ky;b) vs γ(ky;0)): 7.2e-6 rad at ±1e-4, 7.2e-5 rad at ±1e-3;
- winding: wrapped(γ[last]−γ[0]) ≈ 0 for every b — no Chern/winding quantization; any |winding| ≥ π case is **flagged**, never silently smoothed.

Regressions R3/R4 inject a branch-unstable row (π−1e-9) and a raw finite-b drift ≥ π; both must be rejected.

### 1.5 Failure assertions with nonzero exit codes (mandatory)

Every core quantity is now an assertion in `test_multiband_igmc_gate.py` (anchors, b-step, 41² direct/Wilson discrepancy, joint gauge invariance, each negative control, gap scans, branch checks, Liouville) — any violation exits nonzero. `scripts/test_gate_regressions.py` (7 tests, ~10 s) deliberately injects wrong physics (wrong tensor factor, negative indirect gap, branch-unstable row, finite-b drift, plus anchors/pinned-identity/regression numbers) and requires rejection: **7/7 pass, exit 0**; the wrong-tensor-factor test fails by construction (||[H,Sx]|| = 1.20 eV ≥ 1.19 → conservation guard trips).

### 1.6 Independent finite-frequency density-matrix (orbital-Liouville) solver (mandatory)

The ad hoc +0.034304 number is now a **committed, reproducible solver** `orbital_liouville_imag_per_omega`: per-k solve of `[ω I − L] δρ = −[S^x, ρ₀]` on the 16-dimensional Liouville space (ρ₀ = occupied projector), `K(ω) = −⟨Tr(v_x δρ)⟩`, reporting Im K/ω with exact units and quadrature (ω = 0.02 eV, reduced ħ = e = 1). This is a **model** density-matrix response — explicitly **not** a physical orbital-field/Zeeman IGMC.

| N | Im K/ω (Liouville) | direct Dxx/ω | Wilson slope | Liouville vs direct |
| --- | --- | --- | --- | --- |
| 21² | +0.034304 | +0.034294 | +0.033755 | 0.03 % |
| 31² | +0.034304 | +0.034294 | +0.034046 | 0.03 % |

Three independent routes (direct QGT sum, finite-ω density matrix, A-embedded finite-field Wilson) now agree on v2_positive.

### 1.7 Pinned interface content identity (confirmatory)

`scripts/wannier_operator_contract_pinned26.py` = provenance **comment block** + byte-identical upstream module (commit `cdac6418166466b5243e73d8eb291d6dea8d4772`, `papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/wannier_operator_contract.py`). sha256 of the module body = `8343ae0605bfc7748eb978f4ca3c6d612228502e55472a883609404417ded100` = upstream blob sha, asserted in the gate and regression R6. Upstream branch not updated/merged.

## 2. Full-gate summary (regenerated `MULTIBAND_SAME_GAUGE_IGMC_GATE.json`)

- Verdict record: `PASS` with the explicit meaning "gate PASS = allowed to develop a genuine source-authenticated H/S/A file importer; NOT material IGMC, NOT permission for costly DFT/Wannier".
- Two-band anchor: direct 0.0152366470 / Wilson 0.0151827707 — reproduced at printed precision (frozen round3 CSV).
- v2_positive main result unchanged and asserted: direct +0.03429375, Wilson +0.03375450/+0.03404560/+0.03415174 at 21²/31²/41² (1.57 % → 0.72 % → 0.41 %); b-step 1e-4…1e-3 eV stable (0.41 %).
- Joint gauge invariance: direct residual 2.8e-17, Wilson residual 1.1e-12.
- Negative controls (all asserted): NC1 wrong-spin-gauge |Δ| = 0.136; NC2 with-A −0.02536 vs omit-A +0.03429 on v2_spinZ_pos (|Δ| = 0.0597) plus operator-level max|∂H−D_kH| = 0.799 eV·Å; NC3 orbital-pseudospin |Δ| = 0.282; NC4 corrected zero (1.1e-33 / 0.0); NC5 pinned zero baseline (−2.6e-17).
- Commutators: no_soc_corrected 0.0 (all three S_j); no_soc_wrong_factor_bug 1.2 eV; v2_positive ≤ 0.74 eV; pinned ≤ 0.67 eV.
- Fail-closed: gap floor 0.02 eV; indirect gap required positive; Wilson branch jumps/margins/drift asserted; winding reported.

## 3. Test commands (reproduce)

```
cd idea_screen_1/soc_static_1
python3 scripts/test_multiband_igmc_gate.py        # full gate, ~4 min, writes ../MULTIBAND_SAME_GAUGE_IGMC_GATE.json
python3 scripts/test_gate_regressions.py           # fast regressions, ~10 s, exit 0
python3 scripts/audit_source_gate_invariants.py    # researcher-supplied diagnostics, exit 0
python3 -m py_compile scripts/wannier_operator_contract_pinned26.py
```

## 4. Errata (dated)

- 2026-10-23: commit adcd75e `MULTIBAND_SAME_GAUGE_IGMC_GATE.md` NC2 table reversed the with-A/omit-A sign assignment on `v2_spinZ_pos`; source/JSON were correct. Corrected above; JSON re-keyed with explicit labels.
- 2026-10-23: commit adcd75e NC4 claimed a spin-conserving zero using a non-conserving fixture; corrected fixture + commutator assertions replace that claim.
- 2026-10-23: adcd75e gate Markdown stated `-42'm'` row unresolved; unchanged (needs the primary SM PDF — separate source-availability caveat, no new material claim).

## 5. One next step toward authenticated material operators

Build the **fail-closed material-operator importer + matched-model falsification pair**: (i) a schema-validated Wannier90 spinor importer for H(R), spin S(R) and full position connection A(R) + centers from `hr.dat/spn.dat/mmn.dat/amn.dat` of a fresh WAVECAR-retaining AgF₂ B1 run, with hard provenance (input deck SHA, spinor-Wannier settings, gauge checks) and explicit negative controls (spin-gauge mismatch, center double-count, missing operator, unit mismatch) encoded as the same assertion machinery; (ii) before any new DFT, validate the importer on this synthetic gate (it must reproduce +0.0343 on the toy), then re-run this exact gate on the material operators. No new DFT/VASP/PBS/Wannier authorization is requested or implied by this round.
