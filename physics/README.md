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
