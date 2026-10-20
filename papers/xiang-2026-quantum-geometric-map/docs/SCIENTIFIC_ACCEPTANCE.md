# Scientific acceptance (public core, 2026-10-20)

**Overall: PARTIAL.** Software regression PASS verifies only the published code and named numerical protocols, not every claim in Xiang et al., PRB 113 L201406 (2026).

| Branch | Public capability | Scientific gate |
| --- | --- | --- |
| Spin Zeeman geometry | Direct eigenvectors versus analytic two-band response, helicity and U(1) gauge tests | PASS for stated nondegenerate model |
| Four-band spin MNHE | SM Eqs. (16–19) density-matrix currents compared with Eq. (29) | PASS synthetic nondegenerate model; physical driven DC and degeneracy OPEN |
| Tilted Dirac Fig. 2 | Spectrum, local Zeeman geometry, formal Eq. (11), literal SM Eq. (63) branch | PASS for stated model definitions |
| Fig. 2(b) physical ordinate | Literal SI Eq. (11) remains about (2\pi) above the published plateau in independent audits | OPEN |
| Two-color Liouville | Sequential resolvents, both orderings and separately defined equilibrium population derivative | PASS method subset; no general steady-state claim |
| Steady-state static PHE | Collision operator, bath preparation, Ward/contact terms and order of limits | OPEN |
| Orbital OHE | Classical one-band and analytic T=0 Rashba Fermi-ring benchmarks | PASS for named models |
| Orbital PHE (2D Bz) | Literal SM Eq. (56) compared with finite-field semiclassical derivative | PASS model method; printed-line sign difference OPEN |
| Orbital MNHE (strict 2D and 3D) | SBC, Eq. (49) positional shift and literal Eq. (52) comparison | PASS named model subset; factor-two difference OPEN |
| SCBA + ladder | Matched scalar-disorder self-energy and noncrossing Bethe-Salpeter vertex | PASS independent periodic model; not full PHE |
| Four-channel spectral Kubo | (RA+AR-RR-AA), with channel-specific ladders | PASS symmetric finite-eta spectral sector; not full Kubo–Středa |
| Finite-(q,omega) Ward | Matrix resolvent identity with exact finite-difference longitudinal insertion | PASS finite-eta resolvent identity; physical continuity remains OPEN |
| General orbital QMQ/BCQ | General omitted three-band pieces and complete Kubo–Středa physical DC response | OPEN |
| Real Wannier materials | Requires consistent projected magnetic/position operators, Berry connection and validated response protocol | BLOCKED |

This matrix applies only to the code in this GitHub repository. The original v30 research archive includes additional offline programs and results that are not represented by remote CI. In particular, a finite-eta resolvent Ward identity is not a proof of physical charge conservation at finite eta, and a symmetric spectral Kubo sector is not a complete Kubo–Středa calculation.
