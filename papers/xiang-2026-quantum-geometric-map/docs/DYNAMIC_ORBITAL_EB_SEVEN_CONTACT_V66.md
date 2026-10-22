# v66 | Complete projected-position contacts in a finite dynamic mixed E/B model

**Scientific boundary:** **VERIFIED** specified finite-band, finite-frequency response algebra and independently evaluated magnetic derivative; **CONDITIONAL** scalar-gamma bath with prescribed Peierls H(B), projected Pa/Pb including off-diagonal dipoles; **UNRESOLVED** conserving bulk dc magnetotransport and the original Xiang PRB orbital mixed E/B. The original article remains **PARTIAL**. This module is NOT an independently verified implementation of the author's collision operator.

## 1. Complete model and the previously missing electric contact

Given a finite OPEN four-orbital graph, specified real-space Peierls magnetic field B, and **supplied** full Hermitian projected polarization matrices Pa(B), Pb(B), allow [Pa,Pb] ≠ 0. Set hbar=1 and couple a weak electric field via H(B,E) = H_B - E_b Pb(B).

The physical polarization current in this finite projection is **J_a(E,B)=i[H(B,E),Pa(B)]**, including the explicit electric current vertex

\[
J_a^0=i[H_B,P_a], \qquad J_a^E=-i[P_b,P_a].
\]

The latter vanishes for commuting ideal point-orbital positions (v56), but is generically nonzero when the supplied projected positions have off-diagonal dipoles. Omitting it fails the finite-field static equilibrium Ward check. This is a projected-model result; a real material requires all missing remote-band/embedding vertices.

At fixed B and finite gamma>0, frequency omega, chemical potential mu and T>0:

\[
\begin{aligned}
\rho_0&=f_T(H_B-\mu),\\
{\cal L}&=\gamma-i\omega+i\,{\rm ad}_{H_B},\\
\rho_E&={\cal L}^{-1}i[P_b,\rho_0],\\
\sigma_{ab}(\omega,B)&={\rm Tr}(J_a^0\rho_E)+{\rm Tr}(J_a^E\rho_0).
\end{aligned}
\]

Differentiate the **entire** expression in B; retain all seven separately labelled groups:

1. Current magnetic vertex: Tr(i[M,Pa] rhoE), M=dB H;
2. Current magnetic position: Tr(i[H,dB Pa] rhoE);
3. Driving source position: Tr(J L^-1 i[dB Pb,rho0]);
4. Equilibrium preparation: Tr(J L^-1 i[Pb,Df(H)[M]]);
5. Magnetic Liouvillian: Tr(J L^-1{-i[M,rhoE]});
6. Explicit electric current contact derivative: Tr(-i([dB Pb,Pa]+[Pb,dB Pa]) rho0);
7. Electric current contact with equilibrium derivative: Tr(JE Df(H)[M]).

These pieces separately depend on the B-dependent representation and are **not separately observable**. All connection terms of both magnetic current M and magnetic position dB Pa/Pb are indispensable.

## 2. Independent source-free numerical acceptance

Model: explicit four-orbital Peierls graph with independently supplied noncommuting Px/Py matrices, (mu,T,B,q/hbar,gamma,omega)=(0.13,0.37,0.23,-0.84,0.071,0.29), hbar=1.

| Representation | Complex sigma_xy | Complex dB sigma_xy |
| --- | --- | --- |
| symmetric | -0.0002810991671 - 0.0159899938868 i | -0.00640294800294 + 0.00450318857971 i |
| landau_y | same to <3e-13 | same to <3e-13 |
| noncommuting orbital U(B)=exp(iBK) | same to <3e-13 | same to <3e-13 |

The independent **explicit eigenbasis double sum** (not resolvent) agrees with sigma to <1e-12. A **five-point finite-B derivative of the entire independently recomputed finite-frequency response** agrees with the seven-term analytic derivative to at most ~1.1e-12 across three B-dependent frames. In Landau representation the individual seven terms are strongly different from symmetric gauge; deliberately deleting dB Pa, dB Pb or dB JE produces mismatches. Including only a subset is gauge dependent.

A second independent observable-level check avoids the linear-response solver: at omega=0, solve the **finite-E dense Liouville superoperator** \`(gamma+i ad_(H-E Pb)) rho(E)=gamma f(H_B)\`, evaluate **full** \`Tr(i[H-E Pb,Pa] rho(E))\`, then obtain its E/B mixed derivative via a four-corner E±h,B±h finite difference. It agrees with the analytic omega=0 dB sigma in all three frames to <2e-9. This is an independent implementation of *the same scalar-bath model*, **not** an independent conserving microscopic scatterer.

At omega=0, gamma=1e-6 the **complete** static electric response tends to zero: 3.85257e-8 in the chosen model; **omitting JE** leaves 0.05675955, an unphysical leftover for an isolated equilibrium finite graph. The gamma->0 limit does not prove a bulk dc conductor, since all spectra are finite/discrete. The spinless TR-real graph obeys reciprocity sigma_xy(B)=sigma_yx(-B); omega sign reversal gives complex response conjugation.

## 3. Limitations and next v67 gate

- Scalar gamma is a phenomenological bath, **not** a conserving collision kernel. A full response Ward identity and microscopic Keldysh/Kubo–Streda oracle with matching transport/current contacts remain OPEN.
- The full magnetic Wannier position matrix, its field derivative, intercell dipoles and remote-band completion are **inputs**, not constructed for a real material.
- Source Supplemental Material Eq20/31 sign, orbital Eq49/52 factor and omitted three-band channels remain OPEN; this work is not an author-certified source correction.
- No thermodynamic L→∞, omega→0, gamma→0 limit order was established; finite equilibrium response does not imply original nonlinear orbital magnetotransport.

v67 should prioritize **conserving microscopic bath or actual electrode Keldysh comparison** with the identical H,P,M and physically matched observable, not add arbitrary fitted gamma models. Source original Xiang paper **PARTIAL**.

## Reproduction

From the Xiang qgmap paper subproject, install \`python -m pip install -e .\`, run \`python -m unittest discover -s tests -v\`, then \`python scripts/reproduce_v66.py\`. 40 new physics gates plus one independent frozen-numeric test. The reference standalone archive runs v60-v66 previous 178+41 = 219 tests, but remote GitHub Actions separately certifies the integrated repository.
