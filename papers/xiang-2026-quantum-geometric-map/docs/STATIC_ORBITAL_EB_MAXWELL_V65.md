# v65 | Supplied projected position and a static orbital E/B Maxwell oracle

**Scientific status: VERIFIED named finite-equilibrium identities, CONDITIONAL four-site Peierls/projection physics, original Xiang driven orbital mixed E/B UNRESOLVED; source-paper overall PARTIAL.** This is a source-independent method gate, not an erratum or transport calculation.

## Physical model, NOT a transport response

Finite open 4-site Hermitian graph, with explicit positions and Peierls magnetic field B, supplied projected polarization \`P=q diag(y)+D\` with nonzero Hermitian interorbital dipole matrix D. D **must be supplied**: centers alone do not determine full Wannier matrix elements. At T>0 and fixed mu, electric coupling \`H(E,B)=H_B-E P_B\`.

Grand potential \`Omega=-T Tr log(1+exp(-(H-mu)/T))\` yields a **static thermodynamic** mixed derivative:

\[
\partial_E\partial_B\Omega
= -\mathrm{Tr}\big[P_B\,D f_T(H_B)[M_B]\big]
  -\mathrm{Tr}\big[f_T(H_B)\,\partial_B P_B\big],\qquad M_B=\partial_BH_B.
\]

The Fermi derivative \`Df[H](M)\` is evaluated using a spectral divided-difference that remains defined under exact degeneracy. Under a B-dependent magnetic gauge or noncommuting projected-orbital frame \`U(B)\`, *both* \`M_B\` and \`d_B P_B\` acquire derivative connections. The two displayed terms can vary with the chosen basis; their sum cannot. \`U=exp(iBK)\` is only a **basis transformation**; it is not an additional physical magnetic field.

Three **independent full-re-diagonalization** finite-difference routes confirm this response: \`d_B(d_E Omega)\`, \`d_E(d_B Omega)\`, and four-corner \`Omega(E±h,B±h)\`. None uses the Fréchet formula as its finite-difference oracle.

## Canonical ensemble and polarization origin

If \`N(B,mu)\` varies, replacing \`P -> P+c I\` shifts the fixed-mu mixed response by \`-c d_B N\`. This is an **ensemble/origin convention**, not evidence for real origin-dependent magnetoelectric response. For **fixed N**, the physical finite-sample thermodynamic potential is \`F=Omega+mu*N\` with a retuned \`mu(B,E)\`.

The canonical mixed derivative includes

\[
\mu'_B=-N_B/N_\mu,\quad
\partial_E\partial_B F|_N
=-\mathrm{Tr}\big[P(Df(H)[M]+\partial_\mu\rho\,\mu'_B)\big]
-\mathrm{Tr}[\rho\,\partial_B P].
\]

It is independently checked using three finite-difference routes that **numerically refit the chemical potential at each perturbed E/B**, without importing analytic \`mu_B\`. The resulting \`Tr (d_B rho)|_N=0\` Ward identity and invariance under \`P+cI\` are explicit gates.

## Frozen numbers: a dimensionless illustrative 4-orbital graph

At \`mu=.13, T=.37, B=.23, signed q/hbar=-.84\`:

| B-dependent representation | Density/Frechet part, fixed mu | Magnetic position contact | Complete fixed mu | Complete fixed N |
|---|---:|---:|---:|---:|
| Symmetric magnetic gauge | +0.028026623919873 | 0 | **+0.028026623919873** | **+0.038777191760062** |
| Landau magnetic gauge | −0.001646757063044 | +0.029673380982917 | **+0.028026623919873** | **+0.038777191760062** |
| Non-Abelian B-dependent orbital frame | +0.045816710428282 | −0.017790086508409 | **+0.028026623919873** | **+0.038777191760062** |

**Omitting the position contact yields a demonstrably wrong, representation-dependent total.** Shifting polarization by \`2I\` changes fixed-mu total from \`0.028026623920\` to \`0.075326274307\`; fixed-N stays \`0.038777191760\`. For a spinless real/time-reversal-symmetric graph, the static mixed response vanishes at B=0 and is odd under B reversal; a graph with complex loop hoppings has a nonzero B=0 signal.

All numbers are **dimensionless finite-model thermodynamic susceptibilities**, not Hall conductivities or experimental magnetotransport coefficients. The external B response is the defined point-orbital Peierls model and prescribed projected D; realistic B-dependent Wannier basis projection could contain additional terms not supplied here.

## Independent acceptance and next step

40 v65 independent scientific gates exercise magnetic Peierls vertices, full projected dipoles, symmetric/Landau orbital frames, noncommuting B-dependent orbital rotation, **mandatory** \`d_B P\` contact negative control, degenerate Fermi Fréchet, independent Maxwell/Helmholtz finite differences, particle conservation and origin-shift controls, time reversal and fail-closed inputs.

This work **does not** solve the full source manuscript's dynamic orbital magnetoelectric mixed E/B current response, microscopic scattering, Kubo–Středa conserving vertices, magnetic positional shift beyond the supplied finite projection, or its reported omitted three-band and sign disagreements. Original Xiang PRB status stays **PARTIAL**.

**v66 target:** use the validated complete \`P_B,d_B P_B,M_B\` interface within a *dynamic* field-driven response, include all current/source/contact derivatives, and seek a **distinct conserving microscopic response oracle**, not another restatement of the same scalar-relaxation Liouvillian approximation.
