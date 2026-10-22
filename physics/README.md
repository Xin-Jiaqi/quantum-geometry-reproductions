# Physics capabilities | Independently accepted computational primitives

This is the **physics capability registry** for the repository. It is separate in purpose from the original-paper reproduction/history tracked under individual papers. No original paper's possibly incorrect printed formula is made authoritative by this registry. Distinguish **VERIFIED** physical primitives, **CONDITIONAL** protocol-dependent sectors, and **UNRESOLVED** full observables.

| Physics capability | Status | Executable implementation | Independent acceptance |
| --- | --- | --- | --- |
| General nondegenerate N-band metric g_nm, Berry curvature Omega_nm, energy-geometric orbital moment mtilde_n | VERIFIED | [orbital_geometry_v54.py](../papers/xiang-2026-quantum-geometric-map/src/qgmap/orbital_geometry_v54.py) | Spectral operator vs independent eigenprojector finite difference; U1 and constant orbital-unitary invariance; two/three bands |
| Non-Hermitian Zeeman QGT: normal/anomalous metric and curvature; local Dirac-node flux | VERIFIED (nondegenerate model geometry only) | [Cui Flux core and N-band physical-spin solver](../papers/cui-2026-curvature-flux-zeeman-qgt/) | 21 tests: analytic Pauli vs spectral/projector, Wilson winding/flux, spin-orbital basis covariance and degeneracy guard; no finite-frequency/SI/material claim |
| Finite-frequency GHE / reciprocal KME local interband susceptibility and Eq.(16) four-sector expansion | VERIFIED (nondegenerate **interband model** only) | [Cui causal_response.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/causal_response.py) · [scope](../papers/cui-2026-curvature-flux-zeeman-qgt/CAUSAL_RESPONSE_ACCEPTANCE.md) | 39 independent geometry/response tests across Python 3.11 and 3.12; spectral versus N² Liouville, source-frequency sign, static-subtracted reciprocity. No metallic intraband/SI/physical collisions/material KME |
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


## v59 exact 1D finite-band leads and separate periodic bulk charge Kubo

v59 adds exact causal semi-infinite 1D tight-binding chain surface g(E) and dispersive self-energy Σ(E), physical Γ(E), full B derivatives of four-terminal Caroli transmissions and floating Hall voltage, plus independently validated QWZ periodic **charge** Chern number via spectral Berry curvature and Fukui-Wilson links (0,+1,-1,0 for masses -3,-1,+1,+3).

| Capability | Status |
| --- | --- |
| Semi-infinite 1D chain selfenergy, retarded causality and exact finite-band edge | VERIFIED for declared lead Hamiltonian |
| Finite coherent four-probe Hall with B-dependent Σ/Γ and complete derivative | VERIFIED model solver; CONDITIONAL device observable |
| Periodic intrinsic 2D **charge anomalous Hall** QWZ Chern benchmark | VERIFIED independent bulk reference |
| General realistic Wannier electrode, inelastic scattering, full original orbital EB bulk MNHE/PHE | UNRESOLVED |

The optional *B-odd scalar* lead onsite shift is an explicitly prescribed magnetogating sensitivity input, **not** a microscopic time-reversal-compatible orbital lead. Do not equate finite-device voltage and periodic bulk Chern conductivity. [v59 details](../papers/xiang-2026-quantum-geometric-map/docs/INDEPENDENT_LEAD_SELFENERGY_PERIODIC_KUBO_V59.md). Original paper PARTIAL.


## v60: true matrix principal-layer selfenergies and same-H topology

The general n-orbital semi-infinite block-lead retarded Green equation now has an accelerated decimation solver, independently checked by direct **160-layer sparse resolvent inversion**, scalar retarded limits, causal matrix spectral density and lead-orbital unitary covariance. Full **non-diagonal** Caroli Gamma_a and magnetic device-frame contact derivatives (including d_B C and d_B Sigma) are supported for finite 2/4-terminal devices. In a translated four-site graph, omitting d_B C breaks magnetic gauge covariance. **VERIFIED** here means correct operator mathematics for specified B-independent principal layer matrices and retarded eta>0; finite-device Hall voltage is **CONDITIONAL**.

A different same-Hamiltonian test uses identical real-space Fourier hoppings and periodic Bloch QWZ H(k). The **finite open projector's center Chern marker** converges to the **periodic Wilson Chern**, while the finite total marker trace vanishes due to edge cancellation. **VERIFIED geometry/topology only.** This does **NOT** prove a finite Hall voltage equals bulk Kubo charge conductivity or the original source-paper orbital mixed EB.

[Methods and limitations](../papers/xiang-2026-quantum-geometric-map/docs/GENERAL_MATRIX_LEADS_SAME_H_QWZ_V60.md). Original paper science PARTIAL; real Wannier magnetic lead H00(B),H01(B), same-H transport thermodynamic limit and intrinsic orbital EB remain UNRESOLVED.


## v61: exact same-QWZ lattice-to-four-terminal charge Hall bridge

[Implementation](../papers/xiang-2026-quantum-geometric-map/src/qgmap/qwz_device_v61.py) | [28 gates](../papers/xiang-2026-quantum-geometric-map/tests/test_qwz_device_v61.py) | [Physical limits](../papers/xiang-2026-quantum-geometric-map/docs/SAME_QWZ_FOUR_TERMINAL_BULK_EDGE_V61.md).

- **VERIFIED named primitives:** QWZ Fourier hopping identical in periodic Bloch model, finite open device and four microscopic semi-infinite principal-layer leads; coherent Caroli scattering, independent S-matrix unitarity, zero-current floating contacts and phase-gauge covariance.
- **CONDITIONAL charge Hall plateau:** clean T=0, μ=0, η=10^-4, L=3/5/7. m=-1/+1 leads Rxy→+1/-1 h/q², while periodic σxy=+1/-1 q²/h. Relates distinct observables via clean chiral transport and conductivity/resistivity inversion, not via their direct equality. At m=+3, the current is unresolved at measurement threshold and Rxy is **not** reported.
- **UNRESOLVED:** bulk Kubo ↔ general finite-device transport with fully controlled thermodynamic/contact/disorder limits; original full multiband orbital E×B susceptibility. Scientific source audit remains PARTIAL.


## v62: periodic charge Kubo and contact-sensitive local NEGF currents

[Acceptance and explicit negative controls](../papers/xiang-2026-quantum-geometric-map/docs/QWZ_BULK_KUBO_SPATIAL_CHARGE_V62.md). Named analytic Pauli vs eigenstate velocity-Kubo charge Hall σxy(μ,T) and local QWZ strip-lead current continuity/spectral identities: **VERIFIED**. The finite-size contact geometry and measurable Rxy remain **CONDITIONAL**; persistent Hall plateau does not imply invariant local bond current. At metallic μ=1.5, naive reciprocal conductivity is decisively false. Generic Kubo–NEGF physical limit and the original Xiang multiband orbital E×B remain **UNRESOLVED**, original paper status **PARTIAL**.


## v63: non-invasive QWZ voltage-probe limit

The [v63 report](../papers/xiang-2026-quantum-geometric-map/docs/QWZ_NONINVASIVE_PROBE_LIMIT_V63.md) independently scales top/bottom probe amplitudes while holding drive leads fixed. Verified (Sigma,Gammaproptolambda^2) and vanishing local cut-current redistribution as (lambda	o0^+) for the stated finite clean QWZ geometry. The measured weak-probe (R_{xy}=1.00215746,h/q^2) is finite-size/position dependent, not an exact universal value. Generic thermodynamic Kubo–NEGF equality and original multiband orbital (E\times B) remain **UNRESOLVED**; Xiang paper **PARTIAL**.


## v64: exact virtual-limit sensor algebra and finite-coherent caveat

[QWZ two-terminal virtual-probe report](../papers/xiang-2026-quantum-geometric-map/docs/QWZ_VIRTUAL_PROBE_ORDERED_LIMITS_V64.md). New **VERIFIED** model-level primitive: lambda→0+ floating voltage from proposed probe-resolved Caroli spectral weights, using only the two-driver retarded resolvent, independently converged against finite-lambda probe physics with O(lambda²) error. Sensor position can change the limiting Hall reading; size 5–11 at fixed contact width is not a thermodynamic limit. **CONDITIONAL** finite coherent Hall, **UNRESOLVED** original Xiang orbital mixed EB, source audit **PARTIAL**.


## v65: static orbital position-contact Maxwell identity

Named finite thermodynamic primitive [v65](../papers/xiang-2026-quantum-geometric-map/docs/STATIC_ORBITAL_EB_MAXWELL_V65.md): **VERIFIED** fully supplied projected finite-orbital position P with offdiagonal dipoles, complete B-dependent frame derivative d_B P and magnetic vertex d_B H, exact static mixed grand-canonical and fixed-N Helmholtz Maxwell identities, numerical E/B derivatives and number-conserving origin test. Four-site Peierls physics and supplied D remain **CONDITIONAL**. Driven bulk orbital magnetotransport and conserving Kubo-Streda are **UNRESOLVED**; Xiang PRB manuscript still **PARTIAL**.
