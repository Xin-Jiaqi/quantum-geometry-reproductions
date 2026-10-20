# Figure / result acceptance · baseline

Definitions: `PASS` independently demonstrated physical/algorithmic subset; `PARTIAL` some mechanisms verified but numeric reproduction not complete; `BLOCKED` essential data or unresolved physical condition; `NOT STARTED` no actual calculation. `PASS` below will never mean published-pixel matching.

| Target | Frozen equations / parameters | Independent verification | Status at source-map checkpoint |
|---|---|---|---|
| Fig. 1(a) bands | Eq. 15; λ=.6 eVÅ; m=1.1 m_e | analytic `ε±=ℏ²k²/(2m)±λk` and extremum | NOT STARTED |
| Fig. 1(b) Hall C | Eq. 11/13/16; varying μ | geometric occupied-band k sum vs Fermi-surface sum; cutoff + polar convergence | NOT STARTED |
| Fig. 1(c) Z_xy | Eq. 6/15; both bands | analytic Rashba angular texture and gauge rotations | NOT STARTED |
| Fig. 1(d) Z_xx | Eq. 6/15 | odd/even angular symmetries and n↔m antisymmetry | NOT STARTED |
| Fig. 2(a) bands | Eq. 17; β=.25, a=1Å, v=.825 eVÅ, M=.1 eVÅ² | `ε±=±|d|`, Γ gap 0.5 eV | NOT STARTED |
| Fig. 2(b) Hall C | Eq. 11/13/17, same parameters | full BZ T=0 µ-scan, vanishing in gap, periodic mesh convergence | NOT STARTED |
| Fig. 2(c) weighted Q and inset raw Q | Eq. 5/17 | Kubo matrix elements, band gap weighting, bounded-gauge checks | NOT STARTED |
| Fig. 2(d) longitudinal D | Eq. 12/17; frequency printed `10^12 Hz` | finite-frequency AC Kubo imaginary part; ω→0; nonzero inside gap; absolute SI unit audit | NOT STARTED |
| SM Tables I-III | SM §V, Neumann tensor constraints | explicit independent MPG / antiunitary constraints | NOT STARTED |
| SM derivations | §I-II-III-IV | two independent methods, not duplicated code paths | NOT STARTED |
| SM electromagnetic limits | §VI | induced E, frequency, gap assumptions | NOT STARTED |

**Panel gate:** A single Fig. 2(b) matching curve cannot close Fig. 2(a-c-d), let alone §V tensor classification. Published curves are not supplied as numeric data; future digitized values will have provenance and uncertainty.

## Round-1 update after actual execution (2026-10-19)

| Part | Current status | Evidence and boundary |
|---|---|---|
| Fig1(a) | **PASS model** | direct spectrum, band minimum checked |
| Fig1(b) | **PASS model mechanism** | Eq16 numeric polar integral vs analytic 0.06–0.12% difference; no SI baseline |
| Fig1(c-d) | **PARTIAL** | geometric textures and symmetries, near-k=0 plot mask; singular convergence pending |
| Fig2(a) | **PASS model** | Γ gap 0.5 eV and periodic bands |
| Fig2(b) | **PARTIAL** | entire µ-scan shape, exact in-gap zero, selected metal point ≤0.183% 201² vs 301²; peaks not exhaustive |
| Fig2(c) | **PARTIAL** | weighted Q Γ hotspot and separate bare Q plotted; gauge basis test, not numerical author-image comparison |
| Fig2(d) | **PARTIAL** | finite in-gap D and exact-AC Kubo slope, 2π/physical SI normalization not resolved |
| SM Tables I-III | **NOT STARTED** exhaustive | one independent Rashba TR countercheck exposes printed Eq14 ambiguity |
| SM derivations | **PARTIAL** | Liouville/Kubo and Eq13 checked; semiclassical/lesser Green's derivations not independently redone |
| SM VI experimental E | **NOT STARTED** quantitative | geometry/magnetic frequency paper assumptions recorded |

`PASS` means *mechanism/model only*. Full method-level release remains `PARTIAL` and CI/GitHub deployment remains `NOT STARTED`.

## Round-2 SI and antiunitary gate

- **Fig. 1(b): PASS for model absolute SI**: no-fit plateau `-gμB σC_xy = 3.88456977 A/(m T)` for `g=2`, against paper plot ~4. Independently checked main Eq. (16) in SI.
- **Fig. 2(b): PARTIAL absolute SI**: at `μ=1eV`, `401²`, `-gμB σC_xy = -33.49 A/(m T)` vs source curve ~-35; full μ-scan quantitative comparison not digitized.
- **Fig. 2(d): PARTIAL with OPEN PHYSICS DISCREPANCY**: in-gap `+0.017754 A/(m T)`, paper ~`+0.04` (visual estimate), not explained by `101²–401²` meshes or full finite-ω Kubo. No artificial factors added.
- **SM Tables I/II: PARTIAL for three representative MPG constraints**, including separate T parity for C/D; main Eq. (14) common η_T conflicts with main Table I and SM §V.
- SM §VI resolves intended THz *energy* convention (`ℏω≈4meV`) but not the precise author's Fig2(d) numeric code input.

Evidence: `docs/ROUND2_TIME_REVERSAL_AND_UNITS.md` and `results/raw/round2_*`. Full paper status remains `PARTIAL`.

### Round-3 Fig. 2(d) checkpoint (2026-10-19)

- `PASS`: nondegenerate Dirac model displacement kernel Eq.(12) prefactor, direct finite-frequency Liouville, Wilson-loop finite-field polarization derivative, convergence in insulating gap.
- `SOURCE ISSUE`: SM p. 6 Eq. (34) intermediate algebra / `Q` definition carries an apparent ×2 relative to its first line and main Eq.(5).
- `PARTIAL`: original published Fig.2(d) absolute SI ordinate is roughly ×2 compared with independently validated main Eq.(12); precise author data and numerical protocol unavailable. Double-amplitude SM intermediate is diagnostic only.


## Round-4 all-MPG checkpoint

- **PASS for mathematical tensor space:** 122 magnetic point groups, independently reconstructed group closure and 244 9D Neumann spaces from pinned toolkit generators. Longitudinal/Hall/mixed class counts are in `ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`. SM Table I-II selected matrices and source monoclinic y-axis convention explicitly tested.
- **PARTIAL for printed symbol-by-symbol table:** manual equivalence of every SM Table I-II independent parameter symbol (including axis relabeling) not yet completed. Full computer-generated dimension and class data are delivered.
- **Not passed:** external full toolkit Python runtime comparison and author Fig2(d) physical plot ordinate.

## Round-5 Supplemental Tables source audit (2026-10-20)

- SM Table I: 122/122 group cases including zero tensor cases exactly matched at full 9-component subspace level.
- SM Table II: 121/122; group `-42'm'` in p.9 has a missing parameter equality `D_xy=D_yx`, otherwise component support/classification matches.
- SM Table III: 122-group allowed-type classification from round4 still agrees; this round did not manually transcribe every Table III textual token independently.
- Source data, 23 transcription rows, machine comparison, frozen checkpoint in `docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md` and `results/raw/round5_*`. Fig2(d) unchanged, PARTIAL.

### Round 6 independent Fig. 2(d) evidence (2026-10-20)

- A new eigenvector-free Pauli-trace derivation and 2D Gauss-Legendre integration produce `D/ω=0.015234078015630925 Å⁻¹/eV` (N=128²); it agrees with the 301² midpoint result to `2.01e−11`.
- Absolute in-gap D = `0.0177538931944 A/(m·T)` for the explicitly declared cyclic 1 THz, `g=2`; original figure visually appears near `0.04 A/(m·T)` but lacks raw data. No fitted x2 multiplier.
- Model/method: **PASS** on the stated gapped low-frequency protocol. Strict Fig. 2(d) published absolute-ordinate identity: **PARTIAL / SOURCE DATA OPEN**. Metallic regions need their own independent occupation-boundary convergence check.


## Round 7 metallic Fermi-boundary gate (2026-10-20)

`Fig2(b)` **PASS method on paper Eq11/13 metal**, see `ROUND7_METALLIC_FERMI_CONTOUR.md`: direct Fermi-domain integrals with analytic boundary and independent Fermi-line delta integral agree to ≤5.35e−10 Å−1 at seven µ≥0 points; exact hole counterpart follows from Eq17 particle–hole symmetry. `Fig2(d)` **PASS main Eq12 metallic numerical integration**, µ=±1 D=−0.0641784474 A/(m T), but **PARTIAL publisher amplitude** since original plotted data absent and visually approximately ×2. Main text/SM acceptance state for other panels stays unchanged; full status in `FINAL_ACCEPTANCE_MATRIX.md`.
