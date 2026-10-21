# v56 | Independent finite-graph mixed electric-magnetic response

**Goal:** construct source-independent physical response from the verified v55 finite Peierls H(B), not fit disputed Xiang et al. Supplemental Material. Original-paper whole-science audit remains PARTIAL.

## Declared E/B protocol

For a finite OPEN graph with specified positions and hoppings, units hbar=1 and signed q=q_over_hbar, use the actual Peierls H(B), magnetic vertex M=d_B H, electric polarization P_b=q X_b, charge current J_a=i[H(B),P_a] and AC length-gauge perturbation H_E(t)=-E_b(t) P_b. The **specified phenomenological bath** is one gamma>0 relaxing toward rho_eq(B)=f_T(H(B)-mu) at fixed B, mu, T.

First-order driven density:
\[
(\gamma-i\omega+i\,\mathrm{ad}_{H_B})\rho_E=i[P_b,\rho_{eq}],\qquad
\sigma_{ab}=\mathrm{Tr}(J_a\rho_E).
\]

Differentiate the COMPLETE electric response with respect to physical magnetic B, preserving all three terms:
\[
\partial_B\sigma_{ab}=
\underbrace{\mathrm{Tr}(J'_a\rho_E)}_{\text{current contact}}+
\underbrace{\mathrm{Tr}[J_a\mathcal L^{-1}i[P_b,\rho'_{eq}]]}_{\text{equilibrium preparation}}+
\underbrace{\mathrm{Tr}[J_a\mathcal L^{-1}(-i[M,\rho_E])]}_{\text{magnetic Liouvillian}},
\]
where J'_a=i[M,P_a], rho'_eq is the full finite-T Fréchet matrix derivative with respect to H in direction M, and L=gamma-i omega+i ad_H. Each term **depends individually on the B-dependent magnetic gauge**; only their sum is gauge-invariant. A separate double-eigenstate sum validates sigma, and a separate five-point **full finite-B difference** validates the complete derivative without importing its analytic formula.

## Quantitative independent four-site benchmark

Exact v55 open graph, (mu,T,gamma,omega,B,q)=(0.1,0.4,0.08,0.35,0.27,-0.87), hbar=1, transverse Jx/Ey.

| Item | Symmetric gauge | Landau gauge |
| --- | --- | --- |
| Complex sigma_xy(B) | +0.00316501811567 - 0.00441058819912i | same to 5e-14 |
| Complete d_B sigma_xy | -0.00418689839099 - 0.03909337651134i | same to 3e-13 |
| Independent five-point finite B | -0.00418689839102 - 0.03909337651104i | -0.00418689839083 - 0.03909337651090i |
| Current contact | -0.01538258767 - 0.02511886022i | +0.00777232800 - 0.00525650874i |
| Equilibrium Fermi derivative | +0.00933094874 + 0.00277016957i | +0.01277087035 + 0.01440797435i |
| Magnetic Liouvillian | +0.00186474054 - 0.01674468586i | -0.02473009675 - 0.04824484213i |

The current contact is essential: dropping it spoils the independent finite-field derivative and magnetic-gauge covariance. A separate exactly stationary finite-graph equilibrium current is zero even for static finite E and B, so a finite-gamma DC response MUST NOT be mislabeled intrinsic bulk transport.

## Independent physics checks and acceptance

19 model tests: exact v55 Peierls reuse; physical current commutator; independent spectral electric Kubo-like linear response; Liouville residual; two magnetic gauges; independent 5-point B derivative; completeness/contact and gauge tests; finite-T Fréchet f(H) including degenerate spectrum; omega reflection; spinless Onsager symmetry; signed charge and B inversion; equilibrium stationary current zero at finite E/B; gamma→0 discrete-spectrum DC limit; origin shifts; invalid model inputs and fail-closed semantics.

**VERIFIED implementation:** the magnetic current contact, equilibrium matrix Fréchet response, complete magnetic derivative and gauge/Onsager identities *of the declared finite-graph model*.

**CONDITIONAL physics:** actual finite-system AC electric response and B derivative because a scalar gamma relaxation to field-equilibrated f(H) is a *chosen phenomenology*, not a conserving microscopic bath.

**UNRESOLVED:** thermodynamic bulk magnetotransport, real leads and Landauer/Keldysh physics, microscopic scattering/contacts, general Wannier embedding and non-Abelian degeneracies, original-paper disputed equations. No new bulk MNHE/PHE claim.

## Running

From the PRB subproject, call:

    python -m unittest tests.test_finite_graph_eb_v56 -v
    python scripts/reproduce_finite_graph_eb_v56.py
    python scripts/reproduce_finite_graph_eb_v56.py --graph-npz my_finite_open_graph.npz --omega .35

An NPZ input must contain the arrays hoppings (N,N Hermitian) and positions (N,2 real). With no NPZ, the CLI includes the exact 4-site reference. Exact full-precision independent evidence in results/frozen/v56_finite_graph_mixed_EB.json. The local handoff ZIP includes original full source, tests and a graph NPZ example.

**Next:** a true thermodynamic/lead or conserving bath treatment, tested in the same gauge-covariant manner, before accepting any bulk orbital MNHE.
