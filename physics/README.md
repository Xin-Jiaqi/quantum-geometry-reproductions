# Physics capabilities | Independently accepted computational primitives

This is the **physics capability registry** for the repository. It is separate in purpose from the original-paper reproduction/history tracked under individual papers. No original paper's possibly incorrect printed formula is made authoritative by this registry. Distinguish **VERIFIED** physical primitives, **CONDITIONAL** protocol-dependent sectors, and **UNRESOLVED** full observables.

| Physics capability | Status | Executable implementation | Independent acceptance |
| --- | --- | --- | --- |
| General nondegenerate N-band metric g_nm, Berry curvature Omega_nm, energy-geometric orbital moment mtilde_n | VERIFIED | [orbital_geometry_v54.py](../papers/xiang-2026-quantum-geometric-map/src/qgmap/orbital_geometry_v54.py) | Spectral operator vs independent eigenprojector finite difference; U1 and constant orbital-unitary invariance; two/three bands |
| Complete m=l and m!=l decomposition of the local Berry-curvature times orbital-moment **geometric sector** | VERIFIED (that sector only) | Same module, curvature_moment_product | Explicit triple band sum; reference 3-band model and Fermi-level sensitivity |
| Transport coefficient inferred from that geometric sector with a specific RTA | CONDITIONAL | Not a general accepted implementation | Requires explicit physical collision/field protocol and all current terms |
| Original-paper SM Eq49 magnetic positional shift and full orbital MNHE | UNRESOLVED | Remains archived as paper reproduction | Printed 2-versus-1 coefficient difference, omitted three-band paths |
| General orbital E*B magnetotransport including all field contacts, band embedding, degenerate subspaces, Wannier gauge and SI units | UNRESOLVED | No certified general solver yet | Requires a distinct gauge-covariant electromagnetic response oracle |

**Important:** The current implementation lives in the isolated PRB package while these primitives are being validated; the registry marks its capability boundary. When migrated to a shared package, preserve numerical acceptance tests and avoid two divergent implementations.

Quick evaluation of an **external model** (inside the PRB paper subtree): save NPZ arrays hamiltonian (...,N,N) and velocities (2,...,N,N), with v=dH/dk in a **fixed k-independent orthonormal orbital basis**; run the v54 CLI. The output is geometry and a named response **subterm only**, not a measured Hall conductivity.

Original-paper discrepancies (spin Eq20/31, orbital Eq49/52, Fig2 SI) remain documented as evidence, not patched into verified physics. Overall Xiang paper reproducibility remains PARTIAL.


## v55: physical embedding and finite Peierls magnetic vertices

The [v55 physics module](../papers/xiang-2026-quantum-geometric-map/src/qgmap/embedding_magnetic_v55.py) independently verifies:

| Capability | Scope and status |
| --- | --- |
| Covariant velocity d_iH-i[A_i,H] and nondegenerate geometry under smooth non-diagonal orbital W(k) | **VERIFIED**, provided complete physical finite-orbital position A_i |
| Explicit finite OPEN graph straight-path Peierls H(B) and d_B H | **VERIFIED**, with site coordinates, hoppings, signed q/hbar and compatible B units |
| Magnetic-vertex gauge transformation for B-dependent G(B) | **VERIFIED**, includes d_B G contacts; naive omission has model residual 0.522, full transformation ~2.8e-17 |
| Realistic Wannier position inferred only from orbital centers | **CONDITIONAL**, point-orbital center approximation alone is generally insufficient |
| Full orbital E×B conductivity, primitive translation at finite B, degenerate bands | **UNRESOLVED** |

These are physical interfaces, not a certified complete magnetotransport tensor. [Scientific methods](../papers/xiang-2026-quantum-geometric-map/docs/PHYSICAL_EMBEDDING_MAGNETIC_INTERFACE_V55.md). Original paper status PARTIAL.


## v56: Conditional finite-graph mixed E×B response, complete magnetic contacts

Independent v56 physics uses v55 finite OPEN Peierls H(B), P_a=q X_a, J_a=i[H,P_a], and a **declared field-equilibrated scalar gamma bath** in length gauge to calculate the full finite-frequency electric conductivity sigma_ab(B,omega) and its complete magnetic derivative. The magnetic derivative includes current contact d_B J, full equilibrium Fermi-matrix Fréchet derivative, and the magnetic Liouvillian vertex. Individual terms are **magnetic-gauge dependent**, the complete sum alone is invariant; it agrees with a completely independent five-point B difference and a separately coded electric eigenstate spectral sum.

| Capability | Acceptance |
| --- | --- |
| Complete analytic differentiation of finite-model electric response, including current contact, Fermi preparation and Liouvillian terms | **VERIFIED mathematical implementation** |
| Magnetic-gauge equivalence, Onsager, omega conjugation, signed charge reversal and finite-system zero equilibrium current | **VERIFIED model physical checks** |
| Finite open graph mixed E×B prediction with specific gamma bath | **CONDITIONAL**, not intrinsic bulk |
| Complete thermodynamic intrinsic orbital MNHE/PHE and microscopic bath/leads | **UNRESOLVED** |

See [v56 derivation](../papers/xiang-2026-quantum-geometric-map/docs/INDEPENDENT_FINITE_GRAPH_MIXED_EB_V56.md), [code](../papers/xiang-2026-quantum-geometric-map/src/qgmap/finite_graph_eb_v56.py) and [frozen validation](../papers/xiang-2026-quantum-geometric-map/results/frozen/v56_finite_graph_mixed_EB.json). Original Xiang source audit remains PARTIAL.


## v57: independently verified real lead-based Hall device, distinct from v56 bath

Coherent [two- and four-terminal wide-band lead solver](../papers/xiang-2026-quantum-geometric-map/src/qgmap/lead_transport_v57.py) now supplies explicit source/drain leads and self-consistent floating Hall voltage probes. The **finite open graph** uses v55's verified Peierls field, local B-independent Γ_a and exact Caroli transmission. Verified checks: full scattering optical theorem/current conservation, field gauge equivalence of **individual T_ab and d_B T_ab**, Onsager–Casimir reciprocity, two-terminal B-even conductance, 4-terminal B-odd transverse voltage, and independent five-point d_B of the full floating-probe solution.

| Physical scope | Status |
| --- | --- |
| Correct coherent local-wideband NEGF/Caroli and Landauer–Büttiker four-voltage-probe linear algebra, its B derivative and gauge/current identities | **VERIFIED** within explicitly supplied finite model |
| Finite coherent device Hall voltage and B slope for particular Γ,T,site graph | **CONDITIONAL device observable** |
| Scalar-gamma v56 closed-system AC response vs this v57 lead-based DC response | **Different observables; no direct numerical identity claimed** |
| Fixed-density thermodynamic bulk orbital MNHE/PHE, general nonlocal/field-dependent lead self-energy, inelastic response | **UNRESOLVED** |

The 1–4 ring scan is a sparse, changing-geometry *mesoscopic* sequence and **NOT a bulk size-scaling proof**. [Physical derivation](../papers/xiang-2026-quantum-geometric-map/docs/OPEN_LEADS_FOUR_TERMINAL_HALL_V57.md), [frozen output](../papers/xiang-2026-quantum-geometric-map/results/frozen/v57_open_lead_hall_validation.json). Original Xiang paper scientific audit remains PARTIAL.


## v58: fixed-density planar coherent four-probe Hall & numerical-resonance audit

[New physical capability](../papers/xiang-2026-quantum-geometric-map/src/qgmap/planar_hall_scaling_v58.py): genuine Nx×Ny nearest-neighbour square lattice **with fixed site density a^-2, local t, chemical potential and q B a²/hbar per plaquette**, a four-lead *multi-site* wideband self-energy and two explicit contact-scaling choices (per-site gamma or fixed sum per lead). Leads use disjoint centered contacts, full self-consistent Hall voltage probes, and a current-normalized four-probe Hall resistance Rxy=VH/I (h/q²).

**Critical numerical requirement:** 80-point Gaussian energy integration may *fail for narrow coherent resonances*. The default reproducibility script therefore uses **adaptive Fermi-window quadrature split at finite-graph eigenenergies**, with convergence status and independent tighter tolerance. In the 5×5 reference, the unconverged 80-point Hall voltage 0.01452986 and slope -0.06136070 change to **0.01484164 and -0.06470123** after convergence. Fixed Gauss survives only as explicitly **uncertified**.

| Scope | State |
| --- | --- |
| Physically specified multi-contact Caroli calculation, current/gauge/Onsager identities, analytic magnetic derivative and adaptive energy integration | **VERIFIED named finite-coherent-model algorithm** |
| Actual device Hall voltage/resistance for selected contact pattern and Γ scaling | **CONDITIONAL** |
| Bulk thermodynamic orbital MNHE/PHE, fixed electronic carrier density, microscopic disorder/dephasing and real material scattering | **UNRESOLVED** |

Even the 5×5,7×5,7×7,9×5,9×9 fixed-site-density sequence is **not a demonstrated bulk-convergent limit**; the spinless model tests conventional Peierls Hall physics, *not* the disputed source-paper quantum geometry mechanism. [v58 science audit](../papers/xiang-2026-quantum-geometric-map/docs/FIXED_LATTICE_ADAPTIVE_HALL_SCALING_V58.md). Original paper reproduction PARTIAL.
