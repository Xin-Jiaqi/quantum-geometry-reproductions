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
