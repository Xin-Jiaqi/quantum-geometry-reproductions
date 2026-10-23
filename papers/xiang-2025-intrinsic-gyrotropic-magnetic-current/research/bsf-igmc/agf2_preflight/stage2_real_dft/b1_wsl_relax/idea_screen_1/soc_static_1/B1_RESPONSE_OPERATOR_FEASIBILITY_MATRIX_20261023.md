# BSF × IGMC: reusable code versus missing real-material response operators

**2026-10-23 | READ_ONLY_CODE_AUDIT / MATERIAL_RESPONSE_BLOCKED_ON_OPERATORS_AND_EQ14**

Scope: existing public GitHub sources on the PR #38 branch. **No VASP, no material-Wannier construction, no computation authorization.** This is a code/content inspection, not a run of the tests enumerated below.

## 1. Existing tools already in the repository

| Component | Evidence and status | AgF₂ applicability |
| --- | --- | --- |
| Two-band Zeeman geometry | \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/src/igmc_lab/core.py\`, \`pair_kernels\` builds spin–velocity interband matrix products, \`kernel_response\` implements C/D/Kubo for **exactly two nondegenerate bands** | **EXISTS_AND_SCOPE_LIMITED**. Useful as a sign and factor control; not a 12-atom SOC multiband calculation |
| Finite-field Wilson polarization derivative | \`src/igmc_lab/independent_checks.py\`, \`finite_field_wilson_derivative\`; currently a **single-occupied-band two-state Dirac model** | **EXISTS_AND_SCOPE_LIMITED**. Independent oracle; cannot be applied directly to the multi-occupied AgF₂ SOC Hamiltonian |
| Independent full-ω Liouville comparison | \`src/igmc_lab/independent_checks.py\`, \`orbital_liouville_kernel\`; solves **two-level** orbital-density Liouvillian | **EXISTS_AND_SCOPE_LIMITED**; not a material response pipeline or actual orbital-magnetic field coupling |
| Eight-band hopping model plus 4-occupied-band non-Abelian Wilson test | \`research/bsf-igmc/docs/MODEL_AND_VALIDATION.md\` and \`research/bsf-igmc/docs/SCIENTIFIC_ASSESSMENT.md\` | **EXISTS_AND_TOY_ONLY**. The response/Wilson match is documented, but all hoppings and spin SOC couplings are hand-chosen, not Ag/F Wannier-derived; no total-energy double well |
| Magnetic-group tensor algebra | \`src/igmc_lab/tensor_symmetry.py\` supports **two named conventions**, \`channel_resolved_symmetry\` and \`paper_printed_eq14\` | **EXISTS_BUT_PHYSICS_PROTOCOL_AMBIGUOUS**. The alternative is a real unresolved source issue, not a harmless formatting detail |
| Main/SM source map | \`docs/SOURCE_MAP.md\` identifies PRL equations 1–17, 11-page SM, Tables I–III, and flags the Eq. (14) P/T question | **EXISTS_AND_DOCUMENTED**. Primary PDFs are locally archived and not distributed on GitHub; inspecting their exact definitions requires that authorized local source |
| AgF₂ B1 SOC evidence | \`soc_static_1/RUN_RECEIPT.json\`, \`MOMENTS.csv\`, \`BANDS_GAP.json\` | **EXISTS_AS_DFT_DIAGNOSTICS_ONLY**. Site-projected moments and gap are **not** Bloch k-resolved matrices of spin, position, velocity, or orbital magnetic coupling |

## 2. Missing capabilities for quantitative AgF₂ IGMC

| Required artifact | Current status | Falsification/acceptance |
| --- | --- | --- |
| Actual SOC Bloch/Wannier Hamiltonian \(\,H_{mn}(\mathbf k)\,\) for the 12-atom cell | **MISSING** | Must match SOC DFT near and across relevant occupied/empty bands, with tracked disentanglement/gauge and frozen input provenance |
| Full same-gauge physical spin/Zeeman matrix \(\,M^{spin}_{mn,j}(\mathbf k)\,\) | **MISSING** | Verify Pauli/spin convention and \(g\mu_B\) scaling, hermiticity, k-/band-covariance; projected site moments or MAGMOM seed cannot substitute |
| Consistent position/Berry-connection \(A_{mn,i}(\mathbf k)\) or covariant k derivatives | **MISSING** | Handle nonzero intracell positions, Bloch/Wannier gauge, non-Abelian occupied subspace, interband off-diagonal elements and numerical derivatives |
| Material \(\,C_{ij}\,\) and \(\,D_{ij}\,\) conductivity integrals | **MISSING** | Choose unique source/physical protocol, units/ω convention, TR/antiunitary mapping; accept only convergent, non-artifactual values |
| Orbital-B driving / full orbital magnetoelectric response | **MISSING** | LORBMOM is PAW local angular momentum only; orbital field coupling requires its own formulation, not \(\sigma\to L_{\rm PAW}\) replacement |
| Independent finite-B multiband Wilson/response cross-check in *actual* AgF₂ SOC basis | **MISSING** | Must agree in the nonresonant insulating small-field limit with a same-gauge spin-only Kubo slope, after units and factors |
| Partner B2 SOC Hamiltonian / true switching path | **MISSING** | A symmetry-image POSCAR (free to generate) is *not* a proven second structural minimum. Type-II 4%-strained branch distinct |

## 3. Three steps before any new expensive DFT / Wannier work

1. **Resolve the Eq. (14) source/response convention** against the authorized 7-page article and 11-page SM. Main-text Eq. (11) C is a Fermi-surface term, Eq. (12) D has low-ω slope, and Eq. (14) is the disputed tensor prescription. Derive how complex frequency, time reversal, \(J/B\) versus \(J/\dot B\), and real/imag quadratures transform. Log both protocols and pin the physically justified one **without retroactively changing frozen baseline tests**. If unresolved, mark \`BLOCKED_PHYSICAL_CONVENTION\`.
2. **Use only existing toy/model data** to test the specific response component \(\,D_{xx}\,\), low-\(\omega\) extraction, finite-field Wilson cross-check and Eq. (14) sign under \(m'_x\), \(m_z\), \(m'_z\). Distinguish a full complex causal Kubo coefficient from a real coefficient with chosen quadrature convention. Do not claim this establishes a new AgF₂ material signal.
3. **Operator feasibility map**: inspect currently available W90Flow / spinor export code and record which outputs can produce full matched \(H\), spin matrices, center/connection matrices, and k-derivatives. Do not claim \`wannier90_hr.dat\` or Ag site moments suffice. This is a read-only repo audit, not permission to generate or run Wannier.

## 4. A revised, more testable story

The screened AgF₂ SOC moments do **not** defeat \(m'_x\): the time-odd axial operation sends \((M_x,M_y,M_z)\to(-M_x,M_y,M_z)\). This permits the observed weak out-of-plane spin and local orbital moments. Under the previously declared (still convention-conditional) symmetry:

| observable | Type-I \(m_z\) | Type-II \(m'_z\) |
| --- | --- | --- |
| structural \(P_z\) | odd | odd |
| weak axial \(M_z\) / local orbital \(L_z\) | **even** | **odd** |
| bulk insulating displacement \(D_{xx}\) | **odd** | **even** |
| AHE \(z\) pseudovector | even | odd |

**The potentially useful result is a measured, finite IGMC displacement signal reversing independently of a weak out-of-plane magnetic readout, not the fact that a particular magnetic tensor is allowed by symmetry.** SOC B1 confirms a plausible magnetic starting state; it does not prove partner switching or nonzero response.

In an ideal band insulator the conduction IGMC C contribution is the zero-T Fermi-surface term and vanishes in the gap. Hall response may be quantized/nonzero in special insulating phases, but no Hall conductivity or Chern number was computed here. Do not sell C- and D-channel switching as two simultaneously measurable in-gap conduction signals.

## 5. Decision / go-no-go

- **PASS** if the primary-source convention and material operator delivery path are independently audited, and the existing toy response shows a non-artificial, convention-consistent finite D slope.
- **PARTIAL** if Eq. (14) or operator chain remains ambiguous; retain the story as an early hypothesis, do not spend hours building an assumed-response Wannier pipeline.
- **NO-GO / PIVOT** if the physically consistent protocol forbids the selected \(\,D_{xx}\,\), enforces zero for the actual SOC magnetic group, or the source operator interfaces are not realistically implementable. This is not proof that AgF₂ as a material lacks all magnetoelectric effects.

### Historical/ethical controls

The prior B1 SOC screen **exceeded** its 4h authorization by 0.658h; do not treat the disclosure as approval of additional resource use. Any future task must use a local hard wall-clock limiter that survives the loss of chat/SSH monitoring, and obtain fresh explicit user authorization. This read-only study requires none.

**No new DFT job is included in this task or document.**
