# v66 → v67 Research Handoff | 2026-10-22

Repository: Xin-Jiaqi/quantum-geometry-reproductions. Xiang et al. PRB original-paper scientific audit **PARTIAL**.

## New result

v66 builds **finite-frequency scalar-gamma response** for a fully **supplied** noncommuting projected position pair (P_x(B),P_y(B)), and the full magnetic derivatives M=d_BH, d_BPx,d_BPy from v65. H(E,B)=H_B-E_y P_y(B). The model's current MUST include the explicit electric-field vertex Jx(E)=i[H(E,B),Px(B)]=i[H_B,Px]-i E_y[Py,Px]. For noncommuting projected positions the second term is nonzero. This is a finite projected model, not a complete realistic Wannier matrix or the original Xiang source response.

The full analytic d_B sigma_xy at arbitrary finite omega includes **seven terms**: magnetic H current, magnetic Px current, magnetic Py driving source, Fermi equilibrium derivative, magnetic Liouvillian, B derivative of the electric-current contact, and the electric-current contact times the Fermi B derivative.

Finite-model reference: mu=.13, T=.37, B=.23, signed q/hbar=-.84, gamma=.071, omega=.29, four orbitals. Three independent B-dependent magnetic/orbital frame implementations yield sigma_xy=(-.000281099167063-.015989993886819i) and d_B sigma_xy=(-.006402948002937+.004503188579709i). Seven-contact formula vs independent full five-point B derivative residual <1.1e-12; explicit band-pair spectral sigma vs Liouville <1e-12.

**Second independent numerical implementation:** solve full dense steady-state superoperator at finite E, with bath target f(H_B) and current operator including E contact, then four-corner E/B finite differences. At omega=0, computed d_B sigma agrees with analytic v66 scalar-gamma result within ~2e-9 in symmetric, Landau and noncommuting orbital gauges. This is a distinct implementation of the SAME gamma model, not an independently conserving microscopic collision solver.

Crucial negative gate: at omega=0 and gamma=1e-6 the finite projected AC current response with electric current contact is +3.85257e-8, tending to zero; removing it gives +0.05675955 and an incorrect equilibrium limit. Spinless TR-real model obeys sigma_xy(B)=sigma_yx(-B). Without offdiagonal dipoles in ideal point-coordinate P, electric contact vanishes.

## Scientific status and release checks

- **VERIFIED** only named finite-matrix algebra, seven contacts, gauge covariance, spectral/finite-B/finite-E independent finite-model numerical checks.
- **CONDITIONAL** finite open Peierls graph, chosen offdiagonal projected position and **scalar gamma phenomenological bath**.
- **UNRESOLVED** real magnetic Wannier matrix, conserving SCBA/Keldysh/Kubo-Streda microscopic bath, thermodynamic dc order of limits and original Xiang orbital E/B. Original source Eq20/31 signs, orbital Eq49/52 coefficient and omitted 3-band diagrams remain OPEN.
- 40 new v66 independent science tests + 1 frozen-numeric test, v60–v66 independent reference package **219/219 PASS locally**. GitHub integrated QGMap remote CI should separately verify total count and final commit. Frozen GitHub values are compared **numerically** to independent regeneration; do not claim remote file byte equality with reference ZIP unless Git blob checksum matches.

## v67 deep work, scientific priority

1. **Conserving microscopic transport oracle**. Construct explicitly open lead/contact Keldysh with a fully specified H(B), P(B) and physical reservoirs, or conserving SCBA/vertex protocol for the same observable. Test Ward charge conservation, magnetic contact derivatives, (q,omega)->0 order of limits. A second implementation of gamma alone is NOT a conserving oracle.
2. **Complete projected Wannier magnetic embedding**. Supply all magnetic field-dependent position matrix elements and interband/contact terms; do not equate orbital centers to full position.
3. **Return to Xiang source Eq49–52 and printed three-band omissions** with the above fully specified physical model, no forced sign/factor corrections. Only then consider original paper reproduction beyond PARTIAL.

Reference modules: v55 embedding_magnetic; v56 finite_graph_eb (diagonal Pa/Pb); v65 equilibrium_orbital_eb (complete supplied positions, static); **v66 dynamic_orbital_eb (finite omega, electric contact and 7 B-derivative terms)**.
