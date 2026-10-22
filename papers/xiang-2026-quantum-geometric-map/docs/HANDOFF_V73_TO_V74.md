# v73 -> v74 | Two-color finite-q spin Liouville and magnetic preparation

**Date:** 2026-10-22. Original Xiang PRB 113 L201406 (2026) reproduction status **PARTIAL**. Source conventions from original main and SM SHA-pinned v49/v53 audits. v73 is an independent periodic-model method gate, NOT identical to the selected adiabatic tau pole of the publication.

## What was accomplished

The v72 finite-q two-band Bloch scalar-potential/Zeeman problem now has two separately regulated complex external fields wE=omegaE+i etaE and wB=omegaB+i etaB; the output has wSigma=wE+wB. Both orderings are implemented as separate density-matrix amplitudes:

- Magnetic field then electric field: L_Sigma^-1(W R_B(k) - R_B(k+q) W).
- Electric field then magnetic field: L_Sigma^-1(Z(k+q) R_E - R_E Z(k)).
- R_B=L_B^-1[Z,F] and R_E=L_E^-1(W F(k)-F(k+q)W).
- All vertices use exact k→k+q frame overlaps, exact longitudinal finite-q Ward current, and a **declared, not author-certified** transverse midpoint velocity.

Independent implementations use **eigenstate-pair spectral denominators** and **direct 4×4 orbital Sylvester matrix inversions** at every k. Complete and separately ordered responses agree to machine precision; charge Ward residual ~1e-18; nontrivial momentum-dependent gauge covariance passes.

At N=24, q=2π/24, omegaE=.27, etaE=.045, omegaB=.12, etaB=.03, background B0=.11, mu=.32, T=.12, beta=.24, alpha=.13, mass=.32, Wilson=.35, the dimensionless transverse chi(Jx,rho,B)/(-i q) is +.001208260956+.009325386045i (B then E), -.043217183392-.019011725386i (E then B), and **-.042008922436-.009686339341i** (total).

### Crucial static B preparation difference

A *collisionless oscillating Zeeman source* gives R_B=L_B^-1[Z,F] with **zero diagonal band occupation drive**, even as omegaB=0 and etaB→0. Equilibrium pre-preparation in a static field additionally requires f'(epsilon) Z_nn in the Fermi Fréchet derivative. An explicitly **separate opt-in prepared-static-B prescription** adds that piece ONLY for omegaB=0. It converges to the independently differentiated, equilibrium-prepared v72 susceptibility at etaB→0, whereas the collisionless sequence does not.

At etaB=.001/.0001/.00001, discrepancies with independently computed v72 static dB are 2.24949e-4/2.25133e-5/2.25150e-6 with static population, versus 0.0316957/0.0315549/0.0315406 without it. This result is **not** the source paper's dissipative steady-state collision integral or selected electric tau pole. Keeping etaE and etaB independent is crucial.

### Negative controls

Suppressing either two-color time ordering, collapsing the two complex poles, assuming dynamic magnetic response equals reequilibrated static Zeeman derivative, or treating finite nonzero q regulator limit as selected author tau pole is rejected by tests/diagnostics. A finite-q eta scan alone does not establish a homogeneous transport dc limit.

## Acceptance and software

- **VERIFIED** specific toy periodic scalar finite-q/B spin-source Ward; two paths and two regulators; independent eigenpair vs dense orbital Sylvester; static equilibrium magnetic prep vs v72 analytic Fréchet derivative; k frame covariance.
- **CONDITIONAL** momentum finite q, toy two-band, modeled equilibrium B prep, declared transverse current; hbar=1 and dimensionless charge=1.
- **UNRESOLVED** original source SM Eq15 electric selected eta_a pole, Eq20/31 sign, main Eq10/SM Eq62 helicity, full orbital MNHE, complete Wannier intercell positional shift, SM Eq49/52 factor and omitted three-band, original Fig2 normalization and source-body transport current. Source overall **PARTIAL**.
- Source path src/qgmap/two_color_finite_q_spin_v73.py; 45 physical tests tests/test_two_color_finite_q_spin_v73.py; independent frozen regeneration test tests/test_two_color_finite_q_spin_v73_frozen.py; reproduce script scripts/reproduce_v73.py; frozen results/frozen/v73_two_color_finite_q_spin.json; docs/TWO_COLOR_FINITE_Q_SPIN_PREPARATION_V73.md.
- Local initial acceptance **46/46 PASS**. Full integrated GitHub CI and post-merge main must be checked before release.

## v74 research focus

1. Original SM Eq15: *source-exact* selected electric etaE pole coefficient and source-specified E/B orderings in the **homogeneous length gauge**, not a finite q scalar. Keep source sign/velocity conventions and spin-Zeeman Bbar=-g mu_B B explicit. Reproduce source Eq20 or identify exactly which preparation/term blocks it.
2. Prove a valid q→0/omega→0/eta limit and transport magnetization/contact subtraction before connecting any finite-q chi to physical dc bulk PHE. Quantify any noncommutativity.
3. Distinguish dynamical omegaB=0 retarded spin preparation from isothermal static B: physical collision operator needed to pick a steady state, not an invented universal gamma.
4. Keep original SM Eq20/31 and opposite helicity branches explicit. Full orbital MNHE remains a separate Wannier magnetic embedding project including original Eq49/52 factor, three-band omissions, SI normalization.
5. Preserve independent field-derivative, source-authentic and Wannier completeness gates. Never upgrade the overall source paper to COMPLETE until original formula and published numbers close.

**Planning estimates, not physics metrics:** GitHub engineering 90-95%; reusable model-level physics 80%; complete original paper 30-40%. Record actual PR HEAD, final CI/merge SHA, and independent main CI only after verified; keep original source status PARTIAL.


### v73 final release | 2026-10-22, PR, squash and independent main CI certified

- [PR #84](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/84) **squash merged** as `68525556a46d4bd24b778fcaa02eb062f1dbc248`; final tested PR HEAD `ee12ba69d8f495c91ab02d4d50b578778e9b3f6d`.
- Final PR-head [Actions 39779576186](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39779576186): **1053/1053 PASS**, `Ran 1053 tests in 188.041s`, original source-paper status explicitly `PARTIAL`.
- Independent exact squash-merge main [Actions 39780234527](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39780234527): **1053/1053 PASS**, `Ran 1053 tests in 124.236s` on exact commit `68525556a46d4bd24b778fcaa02eb062f1dbc248`; source paper still `PARTIAL`.
- 45 new finite-q two-color spin-source method gates plus one independent frozen-number regeneration. Fresh standalone new v73 **46/46**, previous v72 **38/38**, v71 **43/43 PASS**; local monolithic historical suite exceeded local runtime cap and must **not** be claimed passed locally. Remote full core suite is authoritative.
- **VERIFIED** separately regulated electric/magnetic response insertions and BOTH ordering paths; exact finite-q charge Ward and k-frame covariance; spectral eigenpair vs independent dense Sylvester; separate equilibrium magnetic population channel and its convergence to independent v72 dB.
- **CONDITIONAL** finite-q chosen transverse midpoint current, collisionless input regulators and explicitly opt-in equilibrium static-B preparation. **UNRESOLVED** original author's selected electric 1/eta_E pole, true homogeneous dissipative spin PHE, full orbital magnetic Wannier/MNHE and SM sign/factor/3-band and SI plot issues. Original Xiang paper **PARTIAL**.
- Final documentation-only commits made after verified merge CI modify no v73 science source, tests or numeric frozen evidence. Portable handoff independent snapshot is scientifically tested but not declared byte-identical to every remote Git blob.
