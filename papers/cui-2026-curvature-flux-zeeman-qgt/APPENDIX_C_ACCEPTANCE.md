# Cui et al. Appendix C: conventional spin-current, scientific acceptance (2026-10-22)

**Status: VERIFIED for the specified nondegenerate model-method subset; full paper PARTIAL.**

Source: Cui et al., *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3, main Eqs.(20)–(28) and Appendix C Eqs.(127)–(134). This is an independent reproduction in a synthetic physical-spin × orbital model, not the author's material or a complete measured spin Hall effect.

## The physical current operator and conventions

- Fixed orthonormal, k-independent Bloch basis; H and v_a=d_aH, Hessian K_ab=d_b d_a H, physical spin Pauli sigma_c are all in the same basis. Use dimensionless Pauli operator, with physical spin hbar*sigma/2 understood in O_a,c = {v_a,sigma_c}/4. Natural units e=hbar=1.
- This O is the **conventional spin-current** operator. In a system with spin-orbit torque it is NOT a generally conserved spin current; neither torque-dipole terms nor spin-accumulation/boundary contacts are included.
- Eqs.(127)–(128) intrinsic Kubo: sigma_ab,c^(s) = -2 sum_n f_n Im sum_(m!=n) (O_a,c)nm (v_b)mn/(E_n-E_m)^2 = -sum_n f_n [2 (g_s^A)_ab - (Omega_s^N)_ab].
- For a constant sigma_c, the **exact** band derivative of spin-current expectation is

  partial_b <O_a,c>_n = 2 Re sum_(m!=n) (O_a,c)nm (v_b)mn/(E_n-E_m) + <n|{partial_b partial_a H,sigma_c}/4|n>.

- Therefore, the Fermi-surface RTA term tau E_b int f'_n <O_a,c>_n <v_b>_n is equal to the Fermi-sea form -tau E_b int f_n partial_b <O_a,c>_n on a periodic BZ **only with the complete derivative**. A reduced geometry-only Eq.(133) is exact in this model for strictly momentum-linear H (zero Hessian), **not automatically for a periodic tight-binding model**.

## Three independent controls and frozen numerical evidence

1. Synthetic four-band spin × orbital Hamiltonian (separate Pauli algebra), local k=(0.31,-0.22): direct H+E r finite-field eigenstate derivative = 0.3719631930937383 vs Eq.(127) Kubo = **0.37196319310188974**, for the occupied lowest band spin-current y,z driven along x. Occupation is frozen under the static field; this is a local interband oracle, not a metallic DC experiment.
2. Independent central k finite differences of <O_a,c>_n agree with the spectral derivative including the explicit Hessian in the periodic model; for a strictly linear model the Hessian term vanishes. Full tensor basis covariance and degeneracy/missing-Hessian guards are tested.
3. Midpoint BZ integration (64×64, mu=0.15, T=0.18, tau=1, k∈[-pi,pi)²) for **J_x,z/(tau E_x)**:

| Formula | BZ result (natural units) |
| --- | ---: |
| direct Fermi-surface f' O v | 0.009212602562535185 |
| full Fermi-sea integral including Hessian | 0.009212602562535092 |
| reduced geometric-only term | 0.01833282148550087 |
| missing Hessian correction | -0.00912021892296578 |

The nearly factor-two gap in this **chosen toy model** is a physical omitted-operator contribution, not a universal coefficient or a k-mesh convergence problem. A frozen JSON records these values; the tests prohibit silently replacing the correct full RTA response with geometry-only.

## Reproduction

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
~~~

Previous public package: 39 tests (geometry + causal interband). This PR adds **17** Appendix C tests for 56 total, including E-field rediagonalization, explicit hand-coded spectral sum, spectral-sector recombination, finite k derivative, Fermi surface/sea identity, unitary basis covariance, nontrivial four-band real spin vs trivial two-band Pauli-spin control, and negative controls.

## Explicit exclusions and next actual inputs

This acceptance DOES NOT establish full Eq.(134) for generic Bloch states, non-Abelian Kramers degeneracy, conserved spin current with SOC, scattering/vertex corrections, real-system longitudinal spin transport, material-accurate spin matrices, or SI spin Hall conductivity. The AutoWann H/R results currently lack same-gauge physical spin matrices; no real material result is claimed and no VASP calculation is required to test this module.
