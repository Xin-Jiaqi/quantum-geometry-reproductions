# Physics capabilities | Independently accepted computational primitives

This is the **physics capability registry** for the repository. It is separate in purpose from the original-paper reproduction/history tracked under individual papers. No original paper's possibly incorrect printed formula is made authoritative by this registry. Distinguish **VERIFIED** physical primitives, **CONDITIONAL** protocol-dependent sectors, and **UNRESOLVED** full observables.

| Physics capability | Status | Executable implementation | Independent acceptance |
| --- | --- | --- | --- |
| General nondegenerate N-band metric g_nm, Berry curvature Omega_nm, energy-geometric orbital moment mtilde_n | VERIFIED | [orbital_geometry_v54.py](../papers/xiang-2026-quantum-geometric-map/src/qgmap/orbital_geometry_v54.py) | Spectral operator vs independent eigenprojector finite difference; U1 and constant orbital-unitary invariance; two/three bands |
| Non-Hermitian Zeeman QGT: normal/anomalous metric and curvature; local Dirac-node flux | VERIFIED (nondegenerate model geometry only) | [Cui Flux core and N-band physical-spin solver](../papers/cui-2026-curvature-flux-zeeman-qgt/) | 21 tests: analytic Pauli vs spectral/projector, Wilson winding/flux, spin-orbital basis covariance and degeneracy guard; no finite-frequency/SI/material claim |
| Finite-frequency GHE / reciprocal KME local interband susceptibility and Eq.(16) four-sector expansion | VERIFIED (nondegenerate **interband model** only) | [Cui causal_response.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/causal_response.py) · [scope](../papers/cui-2026-curvature-flux-zeeman-qgt/CAUSAL_RESPONSE_ACCEPTANCE.md) | 39 independent geometry/response tests across Python 3.11 and 3.12; spectral versus N² Liouville, source-frequency sign, static-subtracted reciprocity. No metallic intraband/SI/physical collisions/material KME |
| Conventional physical-spin-current: intrinsic Appendix C Eq.(127)/(128) and complete Hessian-corrected RTA derivative | VERIFIED intrinsic **model subterm**; CONDITIONAL RTA | [Cui spin_current.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/spin_current.py) · [acceptance](../papers/cui-2026-curvature-flux-zeeman-qgt/APPENDIX_C_ACCEPTANCE.md) | 17 new + 39 existing = 56 tests, Python 3.11/3.12 CI; finite-field eigenstate vs Kubo, finite-k band derivative and full Fermi-sea vs surface integral; conventional spin need not be conserved under SOC; no real-material/SI/spin-torque claim |
| Physical-spin vorticity Eq.(27): weighted anomalous spin QGT = half curl of conventional spin current | VERIFIED (nondegenerate fixed-basis **model identity**) | [Cui spin_vorticity.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/spin_vorticity.py) · [scope](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md) | Spectral vs direct eigenstate finite-difference curl and Hessian-aware derivative; mixed Hessian cancels in curl but not symmetric RTA; no macroscopic spin conductivity claim |
| 2D magnetic point-group Jahn character screening for GHE tensor channels | CONDITIONAL (static real-coefficient assumed transformation; **not** validated finite-frequency antiunitary Kubo) | [Cui magnetic_symmetry.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/magnetic_symmetry.py) · [31-generator snapshot](../papers/cui-2026-curvature-flux-zeeman-qgt/data/toolkit_2d_generators.json) | 31 group closures and 124 Reynolds/direct-nullspace agreements; six printed source-Table-I conditional differences retained; cannot infer measured response from invariance alone |
| Microscopic antiunitary finite-frequency GHE/KME selection for TR (1'), C2zT (2'), mxT (m') | VERIFIED (three clean nondegenerate **model** symmetry identities only; not the full 31-group Table I) | [Cui antiunitary_response.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/antiunitary_response.py) · [physical scope](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND15_ANTIUNITARY_MICROSCOPIC.md) | Total 100 tests passing Python 3.11/3.12; exact spinor Hamiltonian/velocity/QGT, full complex frequency spectral Kubo vs orbital Liouville, 2/4-band, momentum/BZ parity, advanced-vs-retarded negative controls. No collisions/material KME/31-group direct microscopic proof |
| All 31 2D MPG generator-to-four-band **local** clean interband response constraints | VERIFIED (31 generic-k nondegenerate **synthetic models**; 62 Jahn candidate ranks; NOT 31 bulk material KME) | [Group microscopic bridge](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/group_microscopic.py) · [scientific acceptance](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND16_ALL31_MICROSCOPIC.md) · [progress](../papers/cui-2026-curvature-flux-zeeman-qgt/REPRODUCTION_PROGRESS.md) | 31 local projective spinor H/V/Kubo covariance residual ≤9.2e-16; 62/62 even/odd in-plane candidate rank matches; 7 Liouville controls; 5 square BZ controls. 118 tests passed Python 3.11 and 3.12 GitHub CI. C3/C6 whole-BZ, Kramers, collisions, real materials and full Table I unresolved |
| Non-Abelian cluster Zeeman QGT + Kramers-safe clean interband response | VERIFIED (exact/near-degenerate **cluster-flat interband model only**) | [nonabelian_response.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/nonabelian_response.py) | 9 tests: exact Kramers doublets, global/subspace unitary covariance, projector-QGT trace/eigenvalue invariants, spectral Kubo vs N² Liouville residual (1.92\times10^{-16}). Near-degenerate clusters are flattened; no intraband metallic/collision KME |
| Triangular-lattice periodic model and first hexagonal Wigner-Seitz BZ | VERIFIED (synthetic model + quadrature symmetry only) | [hexagonal_bz.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/hexagonal_bz.py) | 10 tests; D6-symmetric weighted nodes, exact area (2\sqrt3\pi^2), four-group local covariance < (4.5\times10^{-16}), order-3 tensor projector residual < (2.5\times10^{-18}). Not a converged material response |
| Independent primitive reciprocal-cell BZ integration cross-check for periodic Zeeman response | VERIFIED **numerical integration method only**; physical tensor BZ convergence **UNRESOLVED** | [primitive_bz.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/primitive_bz.py) · [Round19 acceptance](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND19_PRIMITIVE_BZ_CROSSCHECK.md) | Gauss-Legendre vs shifted midpoint primitive torus and D6-symmetrized hexagonal first BZ. 145/145 tests and independent audit on Python 3.11/3.12 (PR #67); group 6 Gauss20 full-even-tensor error vs hex12 2.78e-7 (57% relative). No converged material/SI KME claim |
| k-resolved Zeeman interband hotspot, occupied-transition gap and BZ cancellation diagnostics | VERIFIED (synthetic cluster-interband **method**, physical BZ convergence **UNVERIFIED**) | [bz_hotspots.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/bz_hotspots.py) · [scientific audit](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND20_HOTSPOTS_AND_CANCELLATION.md) · [frozen data](../papers/cui-2026-curvature-flux-zeeman-qgt/results/frozen/round20_khotspot_cancellation.json) | 159/159 tests Python 3.11/3.12; band-pair analytical ±frequency oracle vs prior Kubo; full unprojected BZ, fixed channel and occupancy-threshold negative controls. Group 6 yy even cancels ~98.6% at hexagon order5 (model only); group 6m'm' xy ~21.1%; no converged material KME |
| Exact Fermi occupation-sector partition of clean cluster Kubo (cross-μ / same-side) | VERIFIED (**synthetic interband method only**); BZ response convergence UNRESOLVED | [occupation_partition.py](../papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/occupation_partition.py) · [Round21 acceptance](../papers/cui-2026-curvature-flux-zeeman-qgt/ROUND21_OCCUPATION_PARTITION.md) | PR #74 merged; 167/167 tests pass in Python 3.11/3.12. Pointwise directed pairs reconstruct chi(±z); thermal same-side terms cannot be dropped solely on small Fermi differences. Supplements existing bz_hotspots; no material or SI KME claim |
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

## v66: dynamic projected position and electric-current-contact gate

[Finite AC seven-contact orbital E/B](../papers/xiang-2026-quantum-geometric-map/docs/DYNAMIC_ORBITAL_EB_SEVEN_CONTACT_V66.md). **VERIFIED** specified finite projected P_x,P_y, full B-dependent position and current contacts (including noncommuting electric -i[Pb,Pa]), independent spectral Kubo-like and finite-B derivatives, equilibrium static current cancellation, plus separate finite-E dense superoperator DC oracle. The physical finite-graph scalar gamma relaxation remains **CONDITIONAL** and cannot be described as a conserving microscopic bath or full Kubo-Streda. Real Wannier magnetic embedding and original source paper orbital magnetotransport remain **UNRESOLVED**; Xiang scientific status **PARTIAL**.


## v67: two microscopic semi-infinite leads and differentiated Ward

[Conserving elastic reservoir benchmark](../papers/xiang-2026-quantum-geometric-map/docs/CONSERVING_RESERVOIR_WARD_V67.md). **VERIFIED** exact energy-dependent tight-binding lead surface Green functions; elastic noninteracting on-shell S†S=I and its magnetic derivative; optical theorem, differentiated Caroli/transmission and local physical-site charge continuity with B-dependent site projectors. Finite E_gate/B derivatives include Peierls vertex, B-dependent projected position and contact Sigma/Gamma. This is **CONDITIONAL** coherent finite-device transport with device-only static gate, **NOT** a same-observable v66 or original Xiang orbital Hall Keldysh oracle. Original paper remains **PARTIAL**.

## v68: finite AC micro-lead Keldysh under the v66 projected operator

[Technical report](../papers/xiang-2026-quantum-geometric-map/docs/FINITE_FREQUENCY_SAME_DEVICE_KELDYSH_EB_V68.md). Noninteracting energy-dependent semi-infinite lead Fermi baths yield two-energy lesser Keldysh response of the **same projected device polarization rate operator** as v66; includes electric-current contact -i[Pb,Pa] and magnetic lead/position terms. Fully differentiated Dyson-reservoir charge collision Ward and independent full-field static E/B checks are **VERIFIED finite-model identities**. No phenomenological scalar gamma bath. Device-local field, prescribed projected Wannier dipoles and finite boundary physics remain **CONDITIONAL**; original bulk orbital mixed E/B **UNRESOLVED**, Xiang paper **PARTIAL**.
