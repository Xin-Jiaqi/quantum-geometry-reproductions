# Round36 Source-Scope Audit — Cui arXiv:2604.09725v3

**Target:** Rongjie Cui, Longjun Xiang, Fuming Xu, Jian Wang, *Curvature-Flux
Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum
Geometry*, arXiv:2604.09725v3 (17 Sep 2026), CC BY 4.0.
**Source of truth used here:** the submitted v3 LaTeX (`arXiv:2604.09725v3`,
`dualQGT7.tex`, 1509 lines, downloaded 2026-10-23). Equation labels below are
the v3 labels. The Issue #105 note "v3-versus-v2 equation renumbering" is
honoured by pinning every claim to the v3 text; no v2 numbers are assumed.

This document is the claim-to-test matrix. Every claim of the manuscript
(main text, figures, Table I, appendices) is mapped to an oracle, its
independent implementation, the numeric/analytic target, the error budget and
the current status. Nothing in this file asserts material provenance; the
paper is a model study with no material dataset.

---

## 0. Scope decision

The manuscript is a **two-band 2D Dirac model study** (`H = d(k)·σ`) with:
(i) the fourfold normal/anomalous Zeeman QGT decomposition, (ii) the 2D
Hodge-dual flux/winding representation of local Dirac π₁ topology, (iii) the
four gyrotropic conductivity channels and their reciprocal KME channels, and
(iv) Table I magnetic-point-group selection rules. **It contains no material
data, no .spn archive, no DFPT phonon or EPW scattering tables.** Therefore
Round31–35 (Wannier90 native chain, EPW gate, Fe post-processing) are retained
as *infrastructure* for a future ab-initio extension and are explicitly
excluded from the Round36 acceptance gates.

## 1. Claim inventory and test matrix

Legend: `eq` = equation in v3 text; `impl` = implementation module;
`target` = independent numerical/analytic target; `status` = PASS / PARTIAL /
NOT_IMPLEMENTED / RECORDED-DISCREPANCY. All numeric values are in model
natural units (dimensionless k, energy, σ Pauli matrices).

### 1.1 Main text — geometry (Sections II–III)

| Claim (v3) | Impl | Independent target | Status |
|---|---|---|---|
| Eq.(1) `T^ab_nm = r^a_nm r^b_mn` (conventional QGT) | `core.qgt_eigen/qgt_projector/qgt_analytic` | three routes agree at generic k, band ±; gauge-basis covariance | PASS (frozen Round14–21 tests) |
| Eq.(2) `T^Z,ab_nm = r^a_nm σ^b_mn` (Zeeman QGT) | `zeeman_v3.pair_zeeman_sectors` | eigenstate r,σ matrix elements reproduce analytic two-band identity `T^Z = (d̂×v_a)/(2E) − i v_{a,⊥}/(2nE)` (Appendix A) | PASS (Round14 frozen + Round36 recheck) |
| Eq.(3) fourfold decomposition: `g^{N/A,ab}=(1/2)Re/Im[r^aσ^b+σ^a r^b]`, `Ω^{A/N}=−Re/Im(r×σ)` | `zeeman_v3.pair_zeeman_sectors` | tensor symmetry: g symmetric, Ω vector; normal vs anomalous vanish in Hermitian (σ→r) limit | PASS (Round14 frozen) |
| Eq.(4) `Ω^N_n = (n/2) ∇_k×d̂` | `zeeman_v3.analytic_dhat_sectors` vs eigen | pointwise max abs < 1e-12 at generic k; band-resolved `Σ_n Ω^N_n = 0`, `Σ_n Ω^A_n ≠ 0` | PASS (new Round36 tests) |
| Eq.(5) `g^A,ab_n = −(n/4)[∂_a d̂_b + ∂_b d̂_a]` | ditto | pointwise agreement | PASS (new) |
| Eq.(6) `Ω^A_n = −(1/2)[d̂(∇·d̂) − (d̂·∇)d̂]` | ditto | pointwise agreement | PASS (new) |
| Eq.(7) `g^N,ab_n = (1/4)[(d̂×∂_a d̂)_b + (d̂×∂_b d̂)_a]` | ditto | pointwise agreement | PASS (new) |
| Massive Dirac `Ω^A = −(1/2)(k_x,k_y,2m)/d²`, `∇·Ω^A = −m²/d⁴` | `core.omega_analytic_massive`, `zeeman_v3` | exact at generic k; surface integral `∬_disk ∇·Ω^A = −πR²/(R²+m²)` | PASS (Round14 frozen + Round36) |
| Massless `∇·Ω^A = −(1/2)∇·(k/k²) = C δ²(k)`, `C=−π` | `zeeman_v3.massless_flux_and_winding` | normal line integral `∮ Ω^A·n̂ dl = −π` for several radii/meshes; IR: exclude singular point, flux independent of contour | PASS (Round14 + new) |
| Eq.(10) Gauss index `Q = (1/C)∬ ∇·Ω^A = (1/C)∮ Ω^A·n̂ dl = 1` | `zeeman_v3` | `Q = flux/(−π) = 1` for enclosing contours, both bands | PASS (new) |
| Winding `C_w = (1/2π)∮ ∇θ·dl`, `w=(−k_y,k_x)/k²` | `core.winding_ellipse/winding_field` | `C_w = ±1` for ±chirality, several radii/meshes | PASS (Round14 frozen) |
| Hodge duality `Ω^A = ∗w/2` (planar d_z=0) | `zeeman_v3` | pointwise `(Ω^A_x,Ω^A_y)=(−w_y,w_x)/2`; Eq.(12) finite-difference identities `[∇×Ω^A]_z=(1/2)∇·w`, `∇·Ω^A=−(1/2)[∇×w]_z` | PASS (new) |
| Eq.(13)+ `∮ w·dl = 2πC_w ⇒ ∮ Ω^A·n̂ dl = −πC_w`, `∇·Ω^A = −πC_w δ²` | `zeeman_v3` | massless + massless-massive comparison, chirality sign control | PASS (new) |
| Eq.(10′) `Ω_xy = −(n/2)(∂_x d̂×∂_y d̂)·d̂ = −nm/[2(k²+m²)^{3/2}]`; `m→0: −nπ sgn(m)δ²` | `zeeman_v3.berry_curvature_dirac` | analytic pointwise; Wilson-loop Berry phase `= π(1−m/√(1+m²))` | PASS (Round14 frozen) |
| Local-vs-global distinction: `Ω_xy ≠ ∇·Ω^A` pointwise, integrals agree up to sign as m→0 | `zeeman_v3` | disk integral comparison `∬Ω_xy` vs `∮Ω^A·n̂` at several m | PASS (new) |

### 1.2 Main text — transport (Section IV)

| Claim (v3) | Impl | Independent target | Status |
|---|---|---|---|
| Eq.(exp): `J(ω)=gμ_B Σ_nm∫f_n[iω/ε B×Ω^A − 2iω/ε g^N·B − ω²/ε² B×Ω^N + 2ω²/ε² g^A·B]` | `zeeman_v3.gyrotropic_conductivity` | computed from the primary microscopic expression (Appendix B `supp:eq:sigma`) at small ω; **RECORDED-DISCREPANCY: the O(ω) block of main-text Eq.(exp) carries the opposite sign of the Appendix B derivation (same O(ω²) block); selection rules unaffected** | PASS with recorded discrepancy |
| Eq.(gyroconduc): `σ⁻=(σ_xy−σ_yx)/2`, `σ⁺=(σ_xy+σ_yx)/2` | ditto | tensor antisym/sym split of computed σ_ab | PASS |
| Eq.(rule): `Im σ⁻~ω⟨Ω^A⟩, Im σ⁺~ω⟨g^N⟩, Re σ⁻~ω²⟨Ω^N⟩, Re σ⁺~ω²⟨g^A⟩` | `zeeman_v3` | numeric: extract `Im σ⁻/ω, Im σ⁺/ω, Re σ⁻/ω², Re σ⁺/ω²` at small ω and correlate with the ⟨sector⟩ integrals (sign + slope per the primary derivation; main-text O(ω) sign recorded) | PASS (new) |
| Static-limit handling: dc (ω⁰) terms removed because static B does not induce a current | `zeeman_v3` | `ω→0` limit of σ_ab is regular after declared subtraction of the ω⁰ term; residual quoted | PASS (new) |

### 1.3 Main text — reciprocal KME (Section V)

| Claim (v3) | Impl | Independent target | Status |
|---|---|---|---|
| Eq.(exp1): `M = Σ_nm∫ f_n/(2ε)[iω/ε E×Ω^N + 2iω/ε g^A·E + ω²/ε² E×Ω^A + 2ω²/ε² g^N·E]` | `zeeman_v3.kme_alpha` | computed from the Appendix B final expression (Eq. at tex line 940); **RECORDED-DISCREPANCY: O(ω) block sign opposite between main text Eq.(exp1) and Appendix B** | PASS with recorded discrepancy |
| KME selection rules: `Im α⁻~ω⟨Ω^N⟩, Im α⁺~ω⟨g^A⟩, Re α⁻~ω²⟨Ω^A⟩, Re α⁺~ω²⟨g^N⟩` | ditto | numeric correlation of channel coefficients with ⟨sector⟩ integrals | PASS (new) |
| Onsager reciprocity: GHE↔KME probe same sectors in complementary arrangement | `zeeman_v3` | mapping table σ-channel ↔ α-channel; verify complementary pairing numerically | PASS (new) |

### 1.4 Main text — spin QGT (Section VI) and Appendix C

| Claim (v3) | Impl | Independent target | Status |
|---|---|---|---|
| Eq.(sQGT)+(g_s^{N/A},Ω_s^{N/A}) spin four sectors | `spin_vorticity` (Round14) | spectral pair-weighted spin metric/vorticity; `Ω̄^{A,c}_n=(1/2)∇×v^c_n` | PASS (Round14 frozen) |
| Eq.(eq:spin_current_appendixA)/(B7): spin current in four sectors with τ scaling | `spin_vorticity` spectral RTA | Round14/17 frozen tests; full periodic-lattice RTA Hessian vertex declared NOT included | PASS for spectral part; PARTIAL for generic lattice |
| Appendix C spin-current derivation | `spin_vorticity`, `spin_current` | derivative identity `∂_b v^{a,c}_{nn}` verified | PASS (Round14) |

### 1.5 Table I — 2D magnetic point groups (Section IV, Table I)

| Claim (v3) | Impl | Independent target | Status |
|---|---|---|---|
| 31 2D magnetic point groups; group orders; generator closure | `magnetic_symmetry` + `data/toolkit_2d_generators.json` | orders {1:1,2:2,...,94:12}; orthogonality/closure/XOR tests | PASS (Round14 frozen, 31×4=124 subspaces) |
| Sector transformation rule under each generator class | `zeeman_v3.sector_transformation_rule` | empirical extraction with invariant 2-band probe models: unitary ops (C2, C3, Mx) act on all four sectors with a clean covariant +1 rule (residuals ≤1e-16); **antiunitary ops (C2′, m′) have NO clean ±1 sector rule (residuals 0.09–0.33: the transformation mixes sectors)** — independently confirming the Round-14 "candidate Jahn rules are not dynamical antiunitary theorems" finding | PASS (new) |
| `g^N`, `g^A`, `Ω^N`, `Ω^A` allowed sets (Table I rows) | `zeeman_v3.table1_independent_audit` + `response_symmetry_oracle` | two mechanisms, kept separate: (i) unitary-part Reynolds projection over all 31 groups (valid only for unitary ops; antiunitary groups flagged); (ii) microscopic supp:eq:sigma channels for representative gapped 2-band models of groups 1, 2, m, 2′, m′, mm2, m′m2′, 3 | PASS with recorded rows (see below and §3) |
| Blue-highlighted intersection (all four sectors): `{1,m,2,2′,mm2,m′m2′,3,4,6}` | ditto | **1, 2, m, 2′, mm2: all four channels realized in the 2-band class (ratios vs supp:eq:sigma ≈ 1)** — CONFIRMED. `m′m2′`: Ω^N and g^A realized; Ω^A and g^N are **identically zero in the 2-band class** (the group forces `d_z ≡ 0`, hence `Ω^A_z = g^N_xy ≡ 0` pointwise) — class limitation recorded. `3` (and 4, 6, all C_n, n≥3): Ω^A and g^N realized; Ω^N and g^A are **identically zero in the 2-band class** (C_n covariance forces `d_x+i d_y = z f(k²)` — radial in-plane — whose ε²-weighted BZ integrals vanish; proven + verified on a faithful C3 model) — class limitation recorded | PASS on representative subset; class limitations recorded |
| `g^N` row: `m′` is listed **forbidden** | ditto | the microscopic m′-invariant 2-band model realizes a **nonzero `g^N` channel** (`⟨g^N_xy/ε⟩=0.0137`, Im σ⁺/ω ratio 1.0005) — **RECORDED-DISCREPANCY**: the Table-I `g^N` row omits `m′` yet the antiunitary m′ = m_x·T does not forbid the symmetric transverse channel in the 2-band class (consistent with the Round-14 `m′` conditional mismatch) | RECORDED-DISCREPANCY (new, microscopic) |
| `Ω^N` row "all groups": `m′`, `3`, `4`, `6`, and every C_n-containing group | ditto | within the 2-band class `⟨Ω^N_z⟩` is identically zero for the m′-model and for all C_n (n≥3) models (radial theorem); the paper's "all" row cannot be *realized* by any 2-band model in those groups (multi-band models may realize it — not falsified) | class limitation recorded |
| Gray/T-containing groups `11′, m1′, 21′, mm21′, 31′, 41′, 61′, …` (g^A row, Ω^N "all") | parity argument | every T-containing 2D magnetic group forces either Kramers degeneracy (`11′`) or parity-gapless 2-band models (`m1′, 21′, …`): no gapped 2-band representative exists; the rows are not testable within the 2-band class | class limitation recorded |
| Neumann-principle statement: `3m,4mm,6mm` allow O(ω) Hall (Ω^A); `m′m′2,3m′,4m′m′,4′m′m,6m′m′` allow only symmetric transverse (g^N) | `zeeman_v3.response_symmetry_oracle` (representative subset) | consistent with the unitary-part projection + radial theorem (Ω^A survives C_n+mirror; Ω^N/g^A die by the radial theorem) | PASS on representative subset |

### 1.6 Figures

| Figure | Content | Status |
|---|---|---|
| Fig.1 | fourfold decomposition schematic | structural; the four sectors' existence and vanishing in the Hermitian limit are verified by tests |
| Fig.2(a) | Hodge duality Ω^A vs w | verified numerically (1.1 Hodge rows) |
| Fig.2(b) | GHE/KME channel map | verified numerically (1.2/1.3 rows); main-text O(ω) sign caveat recorded |

### 1.7 Appendices

| Claim | Impl | Status |
|---|---|---|
| Appendix A universal two-band identity, projector-trace form `T^Z,ab_n = i Tr[P_n(v_a·σ)P_ñ σ^b]/Δε`, vector form `−(i/2nE)v_{a,⊥} + (1/2E)d̂×v_a`; d̂-formulas (4)-(7) | `core.qgt_analytic`, `zeeman_v3.analytic_dhat_sectors` | PASS |
| Appendix B Liouville/density-matrix linear response; `supp:eq:sigma`; Eq.(exp)/(exp1) derivations | `zeeman_v3.gyrotropic_conductivity/kme_alpha` | PASS (with §1.2/§1.3 recorded sign discrepancy) |
| Appendix C spin-current derivation | `spin_current`/`spin_vorticity` | PASS (Round14) |

## 2. Numbering, units, gauge and error budget

- **Equation numbering:** v3 labels used verbatim; where the HTML rendering
  renumbers (e.g., unnumbered displayed equations), the TeX label is pinned.
- **Units:** dimensionless model units; σ is the physical (dimensionless)
  Pauli spin, not ħσ/2. Response coefficients are geometric (model natural
  units), not SI; no SI conversion is claimed.
- **Gauge:** all oracle quantities are gauge-invariant (eigenstate or
  projector routes); the response uses the length-gauge interband connection.
- **Error budget:** analytic vs eigenstate residuals target `≤2e-12`
  (generic points), contour integrals `≤2e-3` at 512 points and `≤2e-6` at
  4096 points (Round14 frozen levels), power-law coefficient fits
  `≤1e-6` relative at two ω values, and group-projector tolerances `1e-8`
  (Round14 frozen). The massless distribution statements are verified by
  contour integrals with the singular point excluded (declared IR
  regularization) and by the finite-mass → 0 continuity.

## 3. Recorded discrepancies (not silently corrected)

1. **Main-text Eq.(exp) vs Appendix B (O(ω) sign):** the O(ω) block of the
   gyrotropic current in main-text Eq.(exp) is `+gμ_B[iω/ε(B×Ω^A − 2g^N·B)]`;
   the Appendix B derivation (from `supp:eq:sigma`, tex lines 883–884) gives
   `−gμ_B[iω/ε(B×Ω^A − 2g^N·B)]`. The O(ω²) blocks agree. The selection
   rules Eq.(rule) are sign-independent. The oracle pins signs to the primary
   microscopic expression and quotes the main-text difference explicitly.
2. **Main-text Eq.(exp1) vs Appendix B (O(ω) sign):** same structure for the
   reciprocal KME magnetization (tex lines 936–940 vs 304–308).
3. **Antiunitary sector rule fails (new, Round36):** `sector_transformation_rule`
   shows the four sectors transform under *unitary* spatial ops with a clean
   +1 covariant rule (residuals ≤1e-16) but under *antiunitary* ops (C2′, m′)
   no ±1 rule exists — the transformation mixes the sectors (residuals
   0.09–0.33). Hence any Table-I assignment in antiunitary groups cannot be
   derived from a clean local sector sign rule; it must be decided by the
   microscopic response (done for representatives) or by a future
   finite-frequency antiunitary response theorem. This confirms the Round-14
   "candidate Jahn rules are not dynamical antiunitary theorems" finding with
   explicit numbers.
4. **Table I `g^N` row omits `m′` (new, microscopic contradiction):** the
   m′-invariant 2-band model (m′ = m_x·T, `d=(kx, kx·ky, m+k²)`) realizes a
   nonzero `g^N` channel (`⟨g^N_xy/ε⟩ = 0.0137`, Im σ⁺/ω = −gμ_B⟨g^N_xy/ε⟩
   ratio 1.0005). The source `g^N` row lists `m′` as forbidden. The entries
   `2′` and `m′m2′` (the other two Round-14 conditional groups) are **not**
   contradicted: `2′` realizes all four channels exactly as the blue set says,
   and `m′m2′` realizes Ω^N + g^A with Ω^A/g^N identically zero only because
   the group forces `d_z ≡ 0` in the 2-band class. `m′` is therefore the one
   conditional row with a *microscopic* (not just projector-level) conflict.
5. **2-band-class realization limits (new, recorded as limitations, not source
   errors):** (a) every C_n (n≥3) 2-band model has `d_x+i d_y = z f(k²)` (radial
   in-plane) so `⟨Ω^N_z⟩` and `⟨g^A_xy⟩` vanish identically (θ-integral
   cancellation) — the Ω^N "all" row and the g^A row for 3, 4, 6 (and all
   C_n-containing groups) cannot be realized by 2-band models; (b) `m′m2′`
   forces `d_z ≡ 0` ⇒ Ω^A and g^N vanish identically; (c) every T-containing
   group (`11′, m1′, 21′, mm21′, 31′, 41′, 61′, …`) is either Kramers
   degenerate or parity-gapless in the 2-band class — no gapped representative
   exists. None of these falsify the paper (its classification is for general
   multi-band systems); they delimit what the 2-band `H=d·σ` class can verify.

## 4. Explicit non-targets (kept PARTIAL / excluded)

- No material dataset, .spn archive, DFPT phonon or EPW scattering matrix
  exists for this paper; none is fabricated.
- Round31–35 (Wannier90 native chain, official Fe post-processing,
  EPW/DFPT gates) remain PARTIAL infrastructure for a future ab-initio
  extension; they are not Round36 acceptance criteria.
- The full periodic-lattice spin RTA (explicit Hessian vertex) is not
  claimed (Round14 caveat retained).
- Global BZ convergence of any material band structure: not applicable (no
  material).

## 5. Status

Whole-paper status: **PARTIAL**. All 1.1 geometry rows PASS (Eq. 4–7 verified
pointwise to ≤1e-12; Hodge identities to machine precision; massless
flux/winding and Gauss index PASS). All 1.2–1.3 response rows PASS: the four
GHE channels and the four reciprocal-KME channels are ratio-verified against
the verbatim Appendix-B expressions for every nonzero channel of the
representative groups (GHE ratios 1.0001–1.0038; KME ratios 1.0 to machine
precision; KME Re channels — the ε³-weighted ones — are the cleanest, exactly
1.0 for every model), with the two recorded main-text/appendix O(ω) sign
discrepancies. Table I: the unitary-part projection and the microscopic
response oracle confirm the blue intersection for 1, 2, m, 2′, mm2; the
2-band-class realization limits (C_n≥3 radial theorem, m′m2′ `d_z≡0`,
T-group parity gaplessness) and the **new recorded `m′` g^N contradiction**
(§3 item 4) delimit what the 2-band `H=d·σ` class can verify. Remaining
PARTIAL: a finite-frequency antiunitary response theorem that would decide the
`m′` g^N row (and the other antiunitary-group rows) beyond the 2-band class,
and the future ab-initio extension (Round31–35 infrastructure).
