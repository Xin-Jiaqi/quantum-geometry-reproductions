# IGMC Eq. (14) source-convention audit (Xiang et al., PRL 134, 116301 (2025))

**2026-10-23 ｜ Read-only source audit. Verdict: SOURCE_EQ14 = PARTIAL_RESOLVED** — the paper's self-consistent content (SM Tables I/II transcription + the paper's own model results) fixes the antiunitary convention to **channel-resolved** (η_C = +1, η_D = −1). The main-text Eq. (14) sentence "η_T = +1(−1) for R (RT)", taken literally for both channels, contradicts the paper's own SM tables and its own T-invariant Rashba result; it is documented here as shorthand, not silently "fixed".

## 1. Source availability and hashes

| item | pinned/expected | actual |
| --- | --- | --- |
| Article PDF (7 pp, private archive) | sha256 `6b61f8a1ddd37123053d0de744795f986587d837e2d9706cc5640580d94080c2` | **BLOCKED_SOURCE_UNAVAILABLE** — not found on the reachable Mac filesystem (~/Documents, ~/Downloads, ~/Desktop references tree, hash search) nor on the WSL/Windows machine |
| SM PDF (11 pp, private archive) | sha256 `9efa730976083b892e201659509e75cdc640e88a5065abd190524a99bc82b073` | **BLOCKED_SOURCE_UNAVAILABLE** (same search) |
| Independent copy used | arXiv:2312.04093v2 (journal-ref PRL 134, 116301 (2025), CC BY 4.0), 7 pp main text | sha256 `4b5181f456e1b575db27b1198ec2ae0b3fb63306cc52e62c30ad421c6281f886` — **INDEPENDENT COPY**, not the pinned originals |
| SM | not on arXiv | evidence = repo's visually-checked transcription `src/igmc_lab/paper_table_source.py` (SM Tables I/II, pp. 8–10) at evidence level INDEPENDENTLY_RECOMPUTED (transcription), not CROSS_VERIFIED from raw SM |

All publisher PDFs and SM text: **not uploaded, not committed**; only the short permitted equation definitions and my own derivations are reproduced below.

## 2. Printed main-text definitions (arXiv v2 = PRL version), with anchors

- Eq. (9): `Jα = −gµB Σ fn Znm B cos ωt + ∂t Pα(t)` — conduction term (Zeeman curvature Z, in-phase cos quadrature) + polarization derivative.
- Eq. (10): `Pα = −gµB Σ fn [Nnαβ cos ωt + Fnαβ sin ωt] Bβ` — the sin (out-of-phase) term carries the displacement current; the O(ω²) piece is dropped.
- Eq. (11): `σ^C_αβ = Σ∫_k fn Z^αβ_nm` — conduction IGMC, Fermi-surface property.
- Eq. (12): `σ^D_αβ = Σ∫_k fn (2ω/ε_mn) Q^αβ_nm` — displacement IGMC, Fermi-sea property (ω factor explicit).
- Eq. (13): `σ^C_αβ = −Σ∫_k (∂fn/∂ε_n) v^α_n σ^β_nn` — Fermi-surface identity.
- **Eq. (14)** verbatim (printed): "σ^C/D_αβ = η_T det(R) R_αα′ R_ββ′ σ^C/D_α′β′, where R represents a point group operation, det(R) is responsible for the pseudo nature of σ^C/D_αβ, and **η_T = +1(−1) for R (RT)**." (main text, "Symmetry analysis" paragraph; Eq. (14) is stated for the joint symbol σ^C/D).

**Field-quadrature identification (my derivation).** Writing the AC field B(t) = B0 cos ωt, the current from Eq. (9) is J(t) = −gμB[ΣfZ] B0 cos ωt + ∂tP(t) with ∂tP ∝ ω[·] B0 sin ωt. Under time reversal t → −t: cos ωt is even, sin ωt is odd, J is odd. Therefore the coefficient of the **cos (in-phase) quadrature — σ^C — is time-even**, and the coefficient of the **sin (out-of-phase, ∂tB ∝ sin) quadrature — σ^D — is time-odd**. This is the physical origin of η_C = +1, η_D = −1. (The explicit ω in Eq. (12) is exactly the ∂t factor that makes the D channel time-odd.)

## 3. The convention conflict, with source anchors

**L1 — printed Eq. (14) sentence.** A single η_T = −1 for every antiunitary operation RT, applied to both C and D. Under the group 1' = {1, T} (pure time reversal), this predicts σ^C = 0 and σ^D = 0.

**L2 — SM Tables I/II (transcription).** Row C01 lists groups `1` and `11'` (the latter contains T) with the **general 9-component** C tensor — i.e. T alone does **not** constrain σ^C ⇒ η_C = +1. Row D01 lists `1` and `-1'` (antiunitary inversion) with the general 9-component D tensor — i.e. σ^D is invariant under −1' (which requires η_D = −1, since −1' has det(−I) = −1 and σ' = η_D det(R) R σ Rᵀ = −η_D σ). **The SM tables uniquely fix η_C = +1, η_D = −1.**

**L3 — the paper's own model calculations (recomputed this session with the committed `igmc_lab.core`).**

| model | T-symmetry residual | σ^C | σ^D |
| --- | --- | --- | --- |
| Rashba (Eq. 15), μ = 0 / +0.05 / +0.10 eV | 2.1e-14 (T-invariant) | σ^C_xy = 0.0113 / 0.0113 / 0.0213 (1/Å) **nonzero** | σ^D_xx = 0 at all μ **zero** |
| Dirac (Eq. 17), μ = 0 in gap | 2.07 (T-broken) | σ^C_xy = 1.2e-18 (zero, as the paper states for in-gap μ) | σ^D_xx = +3.04e-4 (eV·1/Å); TR partner (mass_sign=−1) = −3.04e-4 (sign flip) |

The T-invariant Rashba model has nonzero σ^C — this is the paper's own Fig. 1(b). Literal L1 would forbid it under group 1'. The D channel vanishes in the T-invariant Rashba and flips between TR partners — consistent with η_D = −1 (T forbids D).

**L4 — protocol vs SM-transcription projector test (this session).** For all 116 (channel, MPG symbol) pairs transcribed from SM Tables I/II, against the 122-group toolkit generator snapshot:

| protocol | matched pairs | mismatches |
| --- | --- | --- |
| `channel_resolved_symmetry` (η_C=+1, η_D=−1) | **115/116** | 1: group 52 `-42'm'` D-row, predicted dim 1 vs SM dim 2 (see §4) |
| `paper_printed_eq14` (η=−1 both) | 75/116 | 41 — every C-row containing pure T (`11'`, `21'`, `m1'`, `2221'`, `mm21'`, `41'`, `-41'`, …) wrongly predicted **forbidden (dim 0)** while SM says 9/5/4/3/2 components |

**Conclusion.** The literal printed Eq. (14) sentence is not self-consistent with the paper's own SM tables (L2) and its own model results (L3). The internally consistent convention is `channel_resolved_symmetry` (L2 + L4). This is reported as a **convention resolution with both paths' evidence** — not an assertion of a paper error; the pinned SM PDF is unavailable, so the SM leg is capped at INDEPENDENT transcription (BLOCKED for the raw SM text).

## 4. Residual: group 52 `-42'm'` D-row (1/116)

`channel_resolved` predicts D dim 1 ({σ_xy = σ_yx}) for `-42'm'`; the SM transcription row D05 lists pattern `off_xy` (dim 2, σ_xy and σ_yx independent). Both protocols agree with each other here (both dim 1), so this is **not** a C-vs-D convention disagreement; candidate causes are (a) the toolkit snapshot's generator set for group 52, or (b) a SM-transcription nuance for this symbol. Left open, no silent fix; re-verify against a fresh SM page image when the pinned SM becomes available.

## 5. Physical consequences for the AgF₂ B1 m'_x state

For the antiunitary m'_x (R_x = diag(−1,+1,+1), det = −1), channel-resolved protocol gives:

| channel | allowed components under m'_x | note |
| --- | --- | --- |
| C (conduction) | **{σ_xy, σ_xz, σ_yx, σ_zx}** (dim 4) | Hall-like xy allowed (consistent with altermagnet-type xy responses); xx forbidden |
| D (displacement) | **{σ_xx, σ_yy, σ_yz, σ_zy, σ_zz}** (dim 5) | **D_xx (longitudinal displacement) allowed** — the B1 fingerprint |

Protocol robustness: **D_xx is allowed under BOTH protocols** (paper-literal also gives {xx,yy,yz,zy,zz} for D). The **C_xy vs C_xx allowedness differs between protocols** — a falsifiable discriminator for a future response calculation. In the insulator (μ in gap) σ^C ≡ 0 (verified numerically, 1e-18), while the displacement slope remains allowed — exactly the paper's stated Fermi-sea/in-gap picture. **Therefore the B1 "insulating displacement D_xx" story is not forbidden by either protocol; the C-channel Hall IGMC allowedness is protocol-sensitive and must be flagged when quoting SM-based predictions.**

## 6. Evidence levels

- Eq. (9)–(14) transcription + field-quadrature derivation: INDEPENDENTLY_RECOMPUTED (arXiv v2 copy, sha 4b5181f4…).
- SM Tables I/II: INDEPENDENTLY_RECOMPUTED via repo transcription (raw pinned SM PDF BLOCKED_SOURCE_UNAVAILABLE).
- Model tests (Rashba/Dirac): INDEPENDENTLY_RECOMPUTED (committed code, this session).
- Protocol-vs-SM projector test: INDEPENDENTLY_RECOMPUTED (this session).
- Pinned user PDF hashes: NOT_CONFIRMED (not reachable).
