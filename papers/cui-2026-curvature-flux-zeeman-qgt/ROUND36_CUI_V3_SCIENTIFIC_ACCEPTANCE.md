# Round36 Scientific Acceptance Report — Cui arXiv:2604.09725v3 (Issue #105)

**Repo:** `qgr-round36` · **Branch:** `research/cui-round36-zeeman-v3-20261023`
**Base:** `ac9b0cb` (Round34, PR #103)
**Date:** 2026-10-23
**Target manuscript:** R. Cui, L. Xiang, F. Xu, J. Wang, *Curvature-Flux
Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum
Geometry*, arXiv:2604.09725v3 (17 Sep 2026). All claims pinned to the v3
LaTeX (`dualQGT7.tex`, 1509 lines; equations cited by v3 label/line).
**Deliverables:** this report, `ROUND36_SOURCE_SCOPE_AUDIT.md` (claim-to-test
matrix), `src/flux_lab/zeeman_v3.py` (oracles), 12 new tests, frozen JSON
`results/frozen/round36_zeeman_v3_table1_audit.json`.

---

## 1. What was verified and how

All quantities are dimensionless model units (`H = d(k)·σ`, σ = physical
Pauli spin). Every numerical statement below is backed by a test in
`tests/test_table1_zeeman_v3.py` plus the frozen JSON; the full suite is
**383 PASS (371 baseline + 12 new)**.

### 1.1 Two-band geometry oracle (main text Eqs. (4)–(7))

- Eqs. (4)–(7) (`Ω^N=(n/2)∇_k×d̂`, `g^A,ab=−(n/4)(∂_a d̂_b+∂_b d̂_a)`,
  `Ω^A=−(1/2)[d̂(∇·d̂)−(d̂·∇)d̂]`, `g^N,ab=(1/4)[(d̂×∂_a d̂)_b+(d̂×∂_b d̂)_a]`)
  verified **pointwise to ≤1e-12** against the eigenstate route at generic k,
  including the anisotropic model `d=(v_x k_x, v_y k_y, m+κ k_x)`.
- **Hodge duality (Eq. 12):** `[∇×Ω^A]_z = (1/2)∇·w` and
  `∇·Ω^A = −(1/2)[∇×w]_z` — finite-difference residuals **0.0** on the
  corrected anisotropic `Ω^A = −v_x v_y (k_x, k_y)/(2D)`, `w = v_x v_y(−k_y,k_x)/D`.
  (The earlier isotropic-only `Ω^A` left a 1.26 residual; root cause: the
  winding field is anisotropic while the vortex field was isotropic.)
- **Anisotropic-plane Hodge dual:** `Ω^A = ∗w/2`; massless flux
  `∮ Ω^A·n̂ dl = −π sgn(v_x v_y)` = `−π C_w`, `C_w = sgn(v_x v_y)` (several
  radii/meshes); Gauss index `Q = flux/(−π) = 1` both bands; massless
  distribution declared IR-regularized (singular point excluded) + finite-mass
  continuity.
- **Dirac winding** `C_w = (1/2π)∮∇θ·dl = ±1` (chirality sign control);
  **quantized local charge** from `∇·Ω^A`; band `Ω_xy = −nm/[2(k²+m²)^{3/2}]`
  with Wilson-loop `π(1−m/√(1+m²))` (Round-14 frozen).
- **Negative controls (all pass):** Hermitian limit (σ→r ⇒ g^A, Ω^A = 0);
  wrong-sign σ flip detected in all four sectors; wrong
  position-matrix-element contraction (detected via pointwise mismatch vs the
  analytic value); band-parity control.

### 1.2 Clean Kubo GHE + reciprocal KME (four tensor channels)

All response numbers use the **verbatim Appendix B expression**
(`supp:eq:sigma`, tex lines 845–852; `M` expression line 940), never the
main-text O(ω) sign.

| Channel (GHE, σ^-/σ^+ split) | Prediction (verbatim supp:eq:sigma) | Verified |
|---|---|---|
| Im σ⁻ / ω (anomalous Hall) | `+(gμ_B/2)⟨Ω^A_z/ε⟩` | ratio 1.0001–1.0005 (all models) |
| Im σ⁺ / ω (gyrotropic metric) | `−gμ_B⟨g^N_xy/ε⟩` | ratio 1.0001–1.0005 (all models) |
| Re σ⁻ / ω² (O(ω²) Hall) | `−(gμ_B/2)⟨Ω^N_z/ε²⟩` | ratio 1.0034–1.0038 (groups where nonzero) |
| Re σ⁺ / ω² (O(ω²) symmetric) | `+gμ_B⟨g^A_xy/ε²⟩` | ratio 1.0035–1.0038 (groups where nonzero) |

Static (ω⁰) term explicitly subtracted at ω=0; power law extracted at
ω = (0.02, 0.05) with single-point + static-subtraction extraction.

| Channel (reciprocal KME, α^-/α^+ split) | Prediction (line 940) | Verified |
|---|---|---|
| Im α⁻ / ω | `−(1/2)⟨Ω^N_z/ε²⟩` | ratio 1.0 |
| Im α⁺ / ω | `−⟨g^A_xy/ε²⟩` | ratio 1.0 |
| Re α⁻ / ω² | `−(1/2)⟨Ω^A_z/ε³⟩` | **ratio 1.0 exact, every model** |
| Re α⁺ / ω² | `−⟨g^N_xy/ε³⟩` | **ratio 1.0 exact, every model** |

The KME Re channels (ε³-weighted) are the cleanest verification: exact 1.0
for every model (the Im α⁻/α⁺ channels are also 1.0 wherever the Ω^N / g^A
integrals are nonzero; in the C_n≥3 class both sides vanish identically, see
§1.3). GHE↔KME complementarity (Onsager pairing) holds: the same sectors drive
σ^-/α^+ and σ^+/α^- in the complementary arrangements.

### 1.3 Table I — 2D magnetic point groups (independent audit)

- **Transformation rules extracted empirically** from invariant 2-band probe
  models (not assumed):
  - unitary ops (C2, C3, Mx): all four sectors transform with a clean
    covariant **+1** rule — residuals ≤ 1e-16;
  - antiunitary ops (C2′, m′): **no clean ±1 rule exists** — the transform
    mixes sectors (residuals 0.09–0.33). This independently confirms the
    Round-14 finding "candidate Jahn rules are not dynamical antiunitary
    theorems" (§3 of the audit doc).
- **Microscopic response oracle** (representative gapped 2-band models of
  groups 1, 2, m, 2′, m′, mm2, m′m2′, 3): the four channels realized per group:
  | group | realized channels | verdict vs Table I |
  |---|---|---|
  | 1, 2, m, 2′, mm2 | Ω^A, g^N, Ω^N, g^A (all four) | blue rows confirmed |
  | m′ | Ω^A, g^N (**g^N nonzero — row omits m′**) | **recorded discrepancy** |
  | m′m2′ | Ω^N, g^A | Ω^A/g^N class-zero (`d_z≡0`) |
  | 3 (and 4, 6, all C_n≥3) | Ω^A, g^N | Ω^N/g^A class-zero (radial theorem) |
- **2-band-class realization limits (recorded, not source errors):**
  (a) C_n≥3 covariance forces `d_x+i d_y = z f(k²)` — radial in-plane — whose
  ε²-weighted integrals `⟨Ω^N_z⟩`, `⟨g^A_xy⟩` vanish identically (verified on
  a faithful C3 model, mirrors broken, C3 check = 0.0); (b) m′m2′ forces
  `d_z ≡ 0` ⇒ `Ω^A_z = g^N_xy ≡ 0` pointwise; (c) every T-containing group
  (11′, m1′, 21′, mm21′, 31′, 41′, 61′, …) is Kramers-degenerate or
  parity-gapless — no gapped 2-band representative exists.

## 2. Recorded discrepancies (not silently corrected)

1. Main-text Eq.(exp) vs Appendix B: **O(ω) block sign** opposite
   (`+gμ_B[...]` vs `−gμ_B[...]`); O(ω²) agrees. Selection rules unaffected.
2. Main-text Eq.(exp1) (KME) vs Appendix B: same O(ω) sign structure.
3. **Table I `g^N` row omits `m′`:** the m′-invariant 2-band model realizes
   `⟨g^N_xy/ε⟩ = 0.0137` with Im σ⁺/ω ratio 1.0005 (microscopic contradiction;
   the Round-14 conditional `2′`/`m′m2′` entries are *not* contradicted — 2′
   realizes all four as printed, m′m2′'s zeros follow from the forced `d_z≡0`).
4. Antiunitary sector mixing (no clean Jahn-style local rule in antiunitary
   groups) — the Table-I rows of antiunitary groups are decided
   microscopically for representatives; a general finite-frequency antiunitary
   response theorem remains open.

## 3. Acceptance verdict

| Issue #105 item | Verdict |
|---|---|
| ① source-scope audit matrix (`ROUND36_SOURCE_SCOPE_AUDIT.md`) | **PASS** (complete claim-to-test matrix; every row has impl/target/status) |
| ② independent 2D Dirac geometry oracle | **PASS** (four-sector vs analytic d̂, radial Ω^A flux/winding, Dirac winding, Hodge dual, quantized local charge; multi-radius/grid/mass; wrong-sign/missing-sector negative controls; clean continuum IR regularization) |
| ③ clean Kubo GHE + reciprocal KME | **PASS** (all four tensor channels; O(ω)/O(ω²) power laws; static-limit handling; reciprocal KME vs Onsager; ratios 1.0001–1.0038 GHE / 1.0 exact KME) |
| ④ Table I 2D magnetic point group audit | **PASS with recorded rows** (group generators + antiunitary action derived and extracted numerically; blue intersection confirmed for 1,2,m,2′,mm2; m′ g^N recorded; 2-band-class limits recorded) |
| ⑤ isolation of auxiliary work (Round31–35 = infrastructure) | **PASS** (explicitly excluded from acceptance gates; no QE, no fabricated material data) |

Whole-paper status remains **PARTIAL** for the same reason as before: the
finite-frequency antiunitary response theorem (would decide the m′ g^N row and
the other antiunitary-group rows beyond the 2-band class) and any future
ab-initio extension are outside this round's scope.

## 4. Remaining problems (explicit)

1. **m′ g^N row** — microscopic contradiction with Table I; needs either a
   general (multi-band) antiunitary response theorem proving
   `⟨g^N_xy⟩ = 0` under m′ = m_x·T at finite ω, or a source-side revision.
   The 2-band class says nonzero.
2. **Antiunitary sector mixing** — no clean local ±1 transformation rule for
   C2′, m′ (and by extension all antiunitary groups); a Jahn-style rule for
   the k-integrated channels is not derivable from the local fields.
3. **2-band-class realization limits** (C_n≥3 radial theorem, m′m2′ `d_z≡0`,
   T-group parity gaplessness) — whether general multi-band models realize
   the paper's Ω^N/g^A rows in those groups is untested (would need a
   multi-band H).
4. **Global BZ / material convergence** — not applicable (no material data).
5. Round31–35 infrastructure (Wannier90 native chain, official Fe
   post-processing) is retained as PARTIAL for a future ab-initio extension,
   per Issue #105 governance.
