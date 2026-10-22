# v72 | Periodic finite-q spin-Zeeman response, exact charge Ward and physical scope

**Source-grounding:** Xiang et al., PRB 113 L201406 (2026), via the repository's SHA256-verified original PDF/SM audits v49 and v53, and v71 source-observable contact no-go. The authentic source SM uses a **selected electric 1/eta_a pole**, both E/B insertion orderings and published spin/orbital tensors. **This v72 model does not reproduce those complete equations or the source's dc PHE/MNHE. Original paper: PARTIAL.**

## Physical model and the new bulk-periodic observable

Avoid the finite device bounded polarization used in v68–v71. Use a specified **periodic two-band Bloch lattice**

\[
H_0(k)=[\beta\sin k_x+\alpha(2-\cos k_x-\cos k_y)]I
+\sin k_x\sigma_y-\sin k_y\sigma_x+
[m+w(2-\cos k_x-\cos k_y)]\sigma_z.
\]

Apply a **spin Zeeman source** \(H(k,B)=H_0(k)+B\sigma_x\), not orbital Peierls B. Probe a commensurate finite q along y with a scalar potential \(\phi_q\); the declared electric convention is \(E_y=-iq\phi_q\), charge=1, hbar=1. The density vertex between momenta is \(W(k+q,k)=I\) only in a fixed physical orbital basis; under a momentum-dependent frame \(U(k)\) it is \(U(k+q)^\dagger U(k)\). Take transverse current as the *declared* symmetric endpoint velocity vertex \(J_x=\tfrac12[v_x(k+q)+v_x(k)]\). This transverse vertex is **not yet** a fully transport/magnetization-separated original-paper current.

The exact finite-q **longitudinal** lattice vertex is instead
\[
J_y(k,q)=[H(k+q)-H(k)]/q.
\]
The independent band-pair density response is
\[
R_{mn}=\frac{[f(\varepsilon^-_n)-f(\varepsilon^+_m)]\,\langle u_m^+|W|u_n^-\rangle}
{\omega+i\eta+\varepsilon^-_n-\varepsilon^+_m}.
\]
For integrated BZ trace, periodicity with commensurate q enforces \(\langle \mathrm{tr}f(H_-)-\mathrm{tr}f(H_+)\rangle=0\), hence the exact finite-q charge Ward:
\[
(\omega+i\eta)\chi_{\rho\rho}-q\chi_{J_y\rho}=0.
\]

Two **independent implementations**: explicit eigenstate-pair spectral sum versus dense 4×4 orbital Sylvester solve at each k. For the mixed spin response \(\partial_B\chi\), differentiate the Fermi matrix using its nondegenerate/degenerate-safe Fréchet divided differences, and the **full** Liouvillian:
\[
\partial_BR=\mathcal L^{-1}\big[W(\partial_BF_-)-(\partial_BF_+)W
+Z_+R-RZ_-\big],\qquad Z_\pm=U_\pm^\dagger\sigma_xU_\pm.
\]
For a k-independent physical spin Zeeman perturbation only, the current vertices have zero *explicit* B derivative. This omission is **not permitted for orbital Peierls B**.

## Frozen dimensionless numbers and rejection gates

For N=24, q=2π/24, omega=.27, eta=.045, B=.11, mu=.32, T=.12, beta=.24, alpha=.13, m=.32, w=.35:

| finite-q quantity | complex value |
| --- | --- |
| density \(\chi_{\rho\rho}\) | +.006620730104594 − .008283390133627i |
| longitudinal \(\chi_{J_y\rho}\) | +.008251927945343 − .007404839628043i |
| transverse \(\chi_{J_x\rho}/(-iq)\) | +.006895888683488 + .003128753730087i |
| **mixed spin Zeeman** \(\partial_B\chi_{J_x\rho}/(-iq)\) | **−.005924308951739 − .007974278348438i** |

The exact charge and mixed-charge Ward residuals are O(10^-18). Analytic B derivative versus independently rediagonalized **five-point B finite difference** differs by O(10^-13). Spectral versus dense independent Liouville agree at machine precision; all complete responses are invariant under a smooth k-dependent spin frame.

Negative controls: wrongly replacing U(k+q)†U(k) by identity after changing basis creates Ward defect ~8.24e-4, while using an endpoint-averaged velocity in place of exact finite-q longitudinal current creates defect ~1.66e-5. Both are deliberately **failing physical implementations**; they are not treated as alternate valid gauges.

**Scientific classification:** VERIFIED specified periodic finite-q Ward, spin B derivative, frame covariance and numerical methods. CONDITIONAL transverse midpoint vertex, finite eta and q, lattice model. UNRESOLVED original source selected tau pole, dc PHE, two field insertion orders, SM Eq20/31 sign, main Eq10/SM Eq62 helicity, full orbital MNHE Eq49/52 factor and three-band omissions, magnetization transport-current subtraction, real magnetic Wannier position and physical SI Fig2 normalization. **Original article PARTIAL.**

New source, 37 physical tests and 1 independently re-generated frozen-evidence test, reproduction script and frozen JSON. See v72→v73 Handoff. Next: source-identical two-color pole orders, independent bulk response normalization, and no assumption that the finite-q scalar protocol already equals the source dc response.
