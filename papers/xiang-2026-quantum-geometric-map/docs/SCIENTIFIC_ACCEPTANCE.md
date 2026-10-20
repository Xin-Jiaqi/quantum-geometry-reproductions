# Scientific acceptance (public core, 2026-10-20)

**Overall: PARTIAL.** Software regression PASS verifies only the published files and numeric protocol, not every claim in Xiang et al., PRB 113 L201406 (2026).

| Branch | Public capability | Scientific gate |
| --- | --- | --- |
| Spin Zeeman geometry | Direct eigenvectors vs analytic two-band `Z`, correct helicity and U(1) gauge tests | PASS two-band nondegenerate |
| Four-band spin MNHE | `spin_multiband_v5.four_band_audit` computes SM Eqs.(16–19) currents and compares Eq.(29) in a periodic nondegenerate four-band model | PASS synthetic model; physical driven DC and degeneracy OPEN |
| Tilted Dirac Fig.2 | Spectrum, local Z, formal Eq.(11) thermal-step, named printed SM Eq.(63) opposite-sign branch | PASS stated model definitions |
| Fig.2(b) physical ordinate | Literal SI Eq.(11) about `2π` above published plateau in offline independent SI audits | OPEN |
| Two-color Liouville | Generic finite-dimensional sequential resolvents, two orderings and separately defined thermal occupancy derivative | PASS method subset; no uniform-gauge/steady-state claim |
| Steady-state static PHE | Collision operator, bath preparation, Ward/contact terms and order of limits | OPEN |
| Orbital OHE | Classical one-band benchmark; exact T=0 Rashba two-ring and two-helicity analytic formulas | PASS these named models |
| Orbital PHE (2D Bz) | `orbital_phe_v4.orbital_phe_eq56` tests literal SM Eq.(56) against finite-field semiclassical derivative | PASS model method; Eq.(56) sign difference OPEN |\n| Orbital MNHE (strict 2D and 3D) | `orbital_mnhe_v5.orbital_mnhe_2d/3d` tests SBC, Eq.(49) positional shift and literal Eq.(52) | PASS named model subset; Eq.(52) factor-two difference OPEN |\n| General orbital QMQ/BCQ | Full omitted three-band pieces and complete Kubo–Středa/Ward physical DC response are not yet remotely migrated | OPEN/offline only |
| Magnetic symmetry | SM Table I 11 representative tensor nullspaces | PASS for stated generators |
| Real Wannier materials | Need projected spin/orbital magnetic moment and Berry/position operators in a consistent gauge | BLOCKED |

This concise matrix applies to **this GitHub core**. The original 30-round offline research archive contains additional programs/results and separate partial acceptances, but is not represented by the remote CI.
