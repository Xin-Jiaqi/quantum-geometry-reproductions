# AgF₂ × Sliding Ferroelectricity × IGMC
## Research status and Codex handoff (2026-10-23)

**START HERE.** This is the single research entrypoint for local Codex. Repository: [Xin-Jiaqi/quantum-geometry-reproductions](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions). **Active PR: [#106](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/106)**, branch `research/agf2-igmc-real-b1-plan-20261024`. Archived methods/DFT PR: [#38](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38), signed-off head `337eca6` (CLOSED, UNMERGED; retain and do not edit). PR #106 incorporates the earlier work; do not assume it is already on `main`.

**Overall scientific judgment:**
- Xiang et al.'s **model-level IGMC reproduction is substantial but PARTIAL against the original published plots**.
- Sliding-induced spin-Zeeman displacement IGMC is **shown possible in an eight-band toy model**, including connector-specific sign rules; **not a real ferroelectric prediction**.
- AgF₂ **B1 non-SOC relaxation screening and SOC static electronic-state screening pass**; the **real material \(D_{xx}/\omega\) is wholly uncomputed**.
- The synthetic four-band same-gauge H/S/A direct-response, Wilson and Liouville chain and synthetic Wannier90 format import passed tests, **not material-source authentication**.
- The existing PR #106 proposal for a **48-WF / six-hour run is NOT scientifically ready** due to an inconsistent spinor occupied-subspace count and unverified producer/cost assumptions. **No VASP/PBS/Wannier job is authorized.**

This handoff is deliberately **physics first**. Do not spend another iteration beautifying existing scripts or reproducing tests already frozen. The highest-value next step is to determine a physically complete (or honestly controlled truncated) occupied spinor subspace before asking for one real-material calculation.

---

## A. Paper reproduction: what is really established

**Source:** Xiang et al., *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, PRL 134, 116301 (2025), DOI 10.1103/PhysRevLett.134.116301. It concerns a *charge current induced by an AC magnetic field*, not the standard current-induced spin accumulation called the Layer Edelstein effect. The mechanism is based on mixed spin/Zeeman–momentum quantum geometry.

1. **PASS, declared model scope:** Main Eqs. (1)–(13) Q/Z geometries, Rashba analytic conduction channel Eq. (16), lattice Dirac bands and low-frequency displacement channel Eq. (12), plus clean two-band Kubo/Liouville, field-dependent Berry-phase Wilson, independent d-vector/Pauli, finite-frequency wavepacket and lesser-Green checks. Round 8 wavepacket vs lesser Green discrepancy ≤1.4×10⁻¹³ on the 101² two-band test. Clean in-gap Dirac displacement example at f=1 THz, g=2: ≈0.01775389 A/(m·T) under the project's **3D-normalized model convention**. **Not an AgF₂ number.**
2. **PARTIAL / ORIGINAL SOURCE CONFLICTS:** Original Fig. 2(d) absolute amplitude is *not* matched (visual approximately ×2; no original plotting raw data); never rescale an independently verified formula to fit. Main Eq. (14) single antiunitary sign is in tension with channel-resolved real-frequency quadratures. Independent model, transcribed SM tables and 115/116 projector comparisons support \(\eta_C=+1,\eta_D=-1\), but exact publisher SM PDF wasn't reverified in the later audit. One SM group row \(-42'm'\) unresolved. SM Eq. (34) Q factor and Eq. (46) cosine Fourier prefactor were flagged as local source inconsistencies; not proof the author's original plot is erroneous.
3. **122 magnetic point groups:** all classified through independent Neumann and character checks (244 channel subspaces), **243/244** matched to the transcription of the author's source tensor table. This is a substantial group-theory reproduction but is not evidence of a material response.

**Read first, not the entire historical PR timeline:**
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/docs/FINAL_ACCEPTANCE_MATRIX.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/REPORT_ROUND8.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/docs/ROUND8_SEMICLASSICAL_GREEN_CLOSURE.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/IGMC_EQ14_SOURCE_CONVENTION_AUDIT.md`

**Scope guard:** conduction C is a clean zero-T Fermi-surface response and vanishes for chemical potential in a clean insulating gap; displacement D is Fermi-sea/metric-related and can remain finite (AC current ∝ frequency). The physical sign, 2D units and observable phase must be defined, not imported from the 3D paper ordinate.

---

## B. New research idea: symmetry-level result vs actual mechanism

**Hypothesis:** layer sliding and the interlayer magnetic connector may independently control an *insulating displacement magnetic current* and an out-of-plane axial magnetic signal.

In the declared \(\eta_D=-1\) response convention, the candidate connector parity is:

| Quantity upon B1↔B2 | Type I: unitary \(m_z\) | Type II: antiunitary \(m'_z\) |
|---|---|---|
| Structural \(P_z\) | reverses | reverses |
| Displacement IGMC \(D_{xx}\) | **reverses** | **unchanged** |
| Weak out-of-plane axial \(M_z\), local \(L_z\) or AHE \(z\) parity | unchanged | reverses |

This is a **symmetry hypothesis**, not proof the material has a finite response or an accessible electric switching cycle. Literature already established stacking ferroelectric control of magnetic responses; novelty cannot be claimed for combining two published ideas or simply presenting an allowed tensor. Type II is a **different strained magnetic branch** in the AgF₂ source, not freely accessible at zero strain.

**Eight-band constructive proof:** Two layers × two displaced orbitals × physical spin, translation-dependent microscopic two-center hopping without adding a hand-tuned layer-odd onsite term. At ±u=0.19, 31², reduced model units:
- AF-like/unitary connector: \(D_{xx}/\omega=\pm0.003085332739\), indirect gap ≈0.234 eV.
- FM-like/antiunitary connector: \(D_{xx}/\omega=+0.019375728336\) for both u signs, indirect gap ≈0.091 eV.
- Independent four-occupied-band Wilson polarization approaches the geometric result as the link grid increases. AF spin-dependent interlayer tunneling removal reduces D by ~57× (model mechanism probe).
- **Negative result:** model band-energy scan has single sampled translation-identified minimum at u=0.5; ±0.19 are **not** minima. A model dipole proxy is NOT true switchable FE; no double well/barrier has been established.

**Primary evidence:**
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/docs/SCIENTIFIC_ASSESSMENT.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/docs/MODEL_AND_VALIDATION.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/docs/STAGE3_ELECTRONIC_ENERGY_GATE.md`

**Potential publishable story, but NOT yet achieved:** source-grounded **two physically distinct metastable sliding domains** with finite measured/calculated material \(D_{xx}\) switching opposite to the magnetic \(M_z\) parity in Type I, quantified AC-sheet-current size and Faraday electric-field background. For a nonmagnetic bilayer current→spin Layer Edelstein effect, do not conflate the observable or material class with this **magnetic-field→charge-current** IGMC hypothesis.

---

## C. Actual AgF₂ status (do not call it a discovery)

**Source:** C2DB-derived Ag₄F₈ 12-atom bilayer B1 screening candidate, not the original article authors' relaxed POSCAR.

| Stage | Existing evidence | What it does NOT show |
|---|---|---|
| Fixed-cell non-SOC screen | 11 ionic steps, screening force max 0.0782 eV/Å (loose target), Type-I alternating Ag moments ~±0.567 μB, occupation gap 1.3875 eV, nearest interlayer F–F 3.002 Å | Accurate total-energy minimum, true barrier, magnetic ground state |
| SOC `vasp_ncl` static | 64 electronic iterations to EDIFF=10⁻⁶, TOTEN −35.55001581 eV, screened gap 1.3616 eV, Ag \(|S_x|≈0.568\,\mu_B\), local projected \(L_z≈+0.070\,\mu_B\) | SOC easy axis, spinor magnetic-space-group density proof, response magnitude |
| Symmetry | 12/12 source positions obey internal glide \(\{m_x|(1/2,1/2,0)\}\) within 1.28×10⁻⁵ Å; 12 spin/orbital site moment vectors obey conditional \(m'_x\) pairing to printed precision | A distinct B2 ferroelectric domain |
| SOC costs | 8 ranks ×4.658 h ≈37.3 rank-hours; prior authorized 4h cap was **exceeded** and documented | Permission to exceed a later cap |
| Input/output availability | exact source POSCAR/INCAR/KPOINTS, hashes and sanitized receipts available | No saved WAVECAR (LWAVE=.FALSE.), no real Wannier H/S/A, no `D_{xx}`, no B2 partner SCF |

**Important:** antiunitary \(m'_x\) transforms time-odd axial site magnetic moments as \((M_x,M_y,M_z)\to(-M_x,+M_y,+M_z)\); common weak \(M_z\) and Ag projected \(L_z\) are **allowed**. PAW local \(L_z\) is not a physical orbital-Zeeman matrix or orbital IGMC. The unitary B1 internal glide maps **B1 to itself**, not to the opposite FE B2. The latter requires the **source-defined layer sliding / Type-I \(m_z\) connector**, and true local-minima validation.

**Evidence:**
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/IDEA_SCREENING_VERDICT.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/RUN_RECEIPT.json`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/SOC_STATIC_VERDICT.md`
- `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/B1_SOC_MAGNETIC_AND_RESPONSE_FEASIBILITY_AUDIT_20261023.md`

---

## D. Completed method development (reusable; DO NOT rebuild)

From PR #38: synthetic gapped four-band SOC H/S/A joint-gauge direct \(D_{xx}/\omega=+0.0342937475\) and **embedded rank-two Wilson derivative +0.0341517445** on 41² (≈0.414% finite-mesh residual); independent 16-dimensional Liouville referee differs from direct by ≈0.03%; joint constant gauge-invariance and deliberate wrong-operator controls passed. Later tests corrected no-SOC spin-conservation test, global *indirect* gap, Wilson phase winding, Wannier90 U-matrix **column-major** serialization and off-site spin \(S(-R)=S(R)^\dagger\). Synthetic carrier/Wannier90-format parser tests pass. **None is an authenticated ab-initio operator; the formatted spn k-block order cannot be inferred independently from spn alone.**

Read `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/MATERIAL_HSA_EXTERNAL_FORMAT_FIX_20261023.md` and its QA JSON only if working on material H/S/A. Do not regenerate all four-band fixtures, open more parser-fix PRs, or present this synthetic coefficient as the AgF₂ response.

---

## E. Critical current blocker: the proposed 48-WF/6h job is NOT ready

PR #106's `stage3_real_igmc_b1/AGF2_MINIMAL_WANNIER_IGMC_JOB_SPEC_20261024.md` currently proposes **48 WFs and assumes 24 occupied / 24 empty**, while its own evidence estimates **100 valence electrons → 50 occupied spinor bands**. This cannot be the full occupied-band projector. Ag 4d(20)+F 2p(24)+Ag 5s(4) sums to 48 **spatial orbital projectors**, normally **96 spinor Wannier functions** if up/down are both included. The omitted deep valence bands might be legitimately integrated out, but not silently: a genuine low-energy response is a **truncated-sector estimate**, not the full Fermi-sea polarization slope until far-band contributions are bounded. Actual VASP `NBANDS` and per-k band count were not recorded in the public receipt.

The proposed 6h cap is also tight: 72–90 SCF steps extrapolate to 5.24–6.55h **before** Wannier projection from the prior 262s/iteration calibration; changed NBANDS may change this further. Local VASP↔Wannier90 ability to export **actual physical spin `.spn`** remains unverified. `IWANNIER90` in the old job spec is not the documented VASP Wannier90 input tag; check `LWANNIER90`, `LWANNIER90_RUN`, and `LWRITE_SPN` with build evidence. Do not silently correct and launch the job.

The ~26.94 Å vacuum requires a 2D **sheet-normalized** response/units, not a vacuum-dependent bulk coefficient. Future response checks need *denser interpolated k meshes* and signed band/k-resolved cancellation, not just subsampling 12×8 to 8×6.

**Detailed independent review (PR #106):**
`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage3_real_igmc_b1/INDEPENDENT_PHYSICS_REVIEW_BEFORE_COMPUTE_20261023.md`.

---

## F. ONE task local Codex should execute immediately (READ-ONLY, NO VASP)

**Mission: resolve the material occupied-subspace and physical-observable definition, then issue one short GO/NO-GO job revision. Do not continue software development.**

1. **Start at the active branch.** Inspect original frozen B1 SOC OUTCAR/vasprun/EIGENVAL available locally, **READ-ONLY**, if present. Pin actual `NELECT`, `NBANDS`, per-k spinor eigenvalues, gap, and the number of bands in every proposed outer/frozen window. Never reconstruct missing eigenstates from site moments. If original files are unavailable, mark `BLOCKED_RAW_BAND_CENSUS` and use no guessed energy-window band counts.
2. **Make the smallest physical basis decision.** Propose a rigorous full-occupied Wannier manifold with explicit **spinor** projectors and required NBANDS; OR a clearly downfolded low-energy model with a concrete omitted-deep-valence contribution test and no claim to a full Dxx. Explain why it can describe **Fermi-sea Zeeman-induced polarization**, not just reproduce HOMO/LUMO bands. Keep this under a couple of pages.
3. **Verify source producer WITHOUT running it.** Inspect existing installed `vasp_ncl`, `wannier90.x`, relevant compile symbols, saved logs and interface docs for a credible same-source Bloch `mmn/amn/spn/eig` + final `hr/r/U/V` route. If producer is unverified, mark `BLOCKED_SPN_EXPORT`. Do **not** execute any VASP/Wannier binary or create/submit a VASP job.
4. **Propose exactly ONE recalibrated physical compute gate (NOT AUTHORIZED):** work needed, realistic upper wall time based on actual spinor NBANDS, hard-stop independent of agent/SSH, SCF and operator-source prerequisites, one stop condition, k-resolved/2D sheet-normalized signed Dxx/ω output and quantified truncation/cancellation uncertainty. Do not promise B2 or an FE barrier from one B1 point.
5. **Publish only two short deliverables in PR #106:**
   - `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage3_real_igmc_b1/CODEX_BAND_MANIFOLD_READONLY_AUDIT.md`
   - `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage3_real_igmc_b1/CODEX_NEXT_SINGLE_PHYSICS_GATE.md`
   Include exact existing-file references, identified blockers and a single `READY_FOR_SEPARATE_AUTHORIZATION / PARTIAL / NO_GO` decision. Post one concise PR #106 comment. Do not modify PR #38 or past numerical receipts.

**Stop here.** Current user authorization is for code review, text and low-cost read-only analysis only. **No VASP, PBS, Wannier, SOC re-run, B2, NEB, MAE or paid material calculations are authorized**, and no prior approval can be reused. If genuine source data is absent, make the correct BLOCKED decision and stop instead of guessing.

### Future physics decision tree, NOT work in this turn

`reliable B1 full/controlled-subspace Dxx with sign and 2D units` → `physically distinct B2 pair and parity` → `verified two local minima/sliding barrier` → `physical response size, induced Faraday E-background, comparison to prior BSF/AHE/LEE`. Conversely `Dxx ≈0` or large uncontrollable deep-band cancellations → analyze and demote AgF₂ candidate; do not endlessly optimize the same toy.

**What a successful paper would have to show:** a *finite, gauge-stable, physically meaningful real-material response* whose switching is tied to a **real bistable structural pathway** and whose magnetic-connector parity supplies a genuinely independent experimental signature. None of this is established yet.
