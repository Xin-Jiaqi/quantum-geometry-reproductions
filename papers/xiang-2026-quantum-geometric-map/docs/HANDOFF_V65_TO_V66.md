# v65 -> v66 | Exact finite-static orbital-position Maxwell handoff

Date: 2026-10-22. Repo: Xin-Jiaqi/quantum-geometry-reproductions.

## New physical result (v65)

Source-independent, finite-open four-orbital graph with explicitly supplied **off-diagonal projected position dipoles**. Magnetic B enters through a specified v55 finite Peierls model and the projected polarization is P=q diag(position_y)+D with the supplied Hermitian D. B-dependent symmetric↔Landau and noncommuting orbital basis transformations must include both M=d_B H and d_B P. The fixed-mu equilibrium thermodynamic mixed Maxwell coefficient is -Tr[P L_f(H;M)]-Tr[f(H)d_B P]. This is STATIC, not a driven transport current. Fermi Fréchet is stable under exact degeneracies.

At mu=.13,T=.37,B=.23,q/hbar=-.84, symmetric, Landau and noncommuting orbital representations give the same full fixed-mu value +0.028026623919873. Individual fixed-mu first/contact terms are (sym +.028026623920, 0), (Landau -.001646757063,+.029673380983), (orbital +.045816710428,-.017790086508). Omitting d_B P fails the magnetic gauge/basis test. Three independent full re-diagonalization finite-difference Maxwell routes validate the sum.

Under fixed N, an analytically derived mu_B=-N_B/N_mu enters, and a physically distinct canonical Helmholtz mixed derivative is +0.038777191760062. Three independently re-fitted chemical-potential finite-difference routes validate it. A constant shift P->P+2I changes the grand-fixed-mu derivative to +0.075326274307 but leaves fixed-N exactly unchanged. Fixed-N electron-number derivative residual ~2e-17. A TR-real spinless graph gives zero B=0 and B-odd response. These are essential ensemble, contact and symmetry acceptance gates.

## Boundaries and unresolved science

VERIFIED only for fully specified finite-T finite-graph equilibrium mixed Maxwell and supplied projected dipoles. The Peierls graph and dipole matrix D are CONDITIONAL model choices, NOT derived from a real Wannier90 material with all B-dependent matrix elements. This is NOT a Kubo-Streda driven dc orbital Hall or the original Xiang PRB mixed E/B response. Original manuscript sign, factor-2 and three-band issues remain unresolved; overall Xiang source science PARTIAL.

This must not be rebranded "full orbital response", "complete Wannier embedding", "bulk thermodynamic DC" or "author method reproduced". Finite static grand-potential derivative does NOT substitute for the original driven current response.

## Files and tests

- src/qgmap/equilibrium_orbital_eb_v65.py
- tests/test_equilibrium_orbital_eb_v65.py (40 independent named gates)
- scripts/reproduce_v65.py and results/frozen/v65_static_equilibrium_orbital_eb.json
- docs/STATIC_ORBITAL_EB_MAXWELL_V65.md, SCIENTIFIC_ACCEPTANCE.md
- v60–v64 standalone baseline 137 + new v65 40 = 177 scientific local tests, distinct from full GitHub CI.
- At release, record exact PR and merged-main GitHub Actions status/count, commit SHA and scope; do not infer remote acceptance from local tests alone.

## Recommended v66 project

Use **the same fully supplied** H_B,M_B,P_B,d_B P_B in the existing v56 **dynamic** finite AC Liouvillian EB response, supplement BOTH the explicit current vertex derivative and the B derivative of electric source from d_B P. Test full analytic derivative independently against five-point finite B differentiating the complete AC response in two gauges/noncommuting orbitals. This is still conditional on gamma; the **next additional gate** must be an independent microscopic conserving bath/current oracle or well-defined thermodynamic scattering construction. Only after such an oracle could paper-level orbital transport be reconsidered.
