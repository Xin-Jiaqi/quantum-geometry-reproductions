# v68 | Same projected-device current in finite-frequency reservoir Keldysh

**Acceptance: VERIFIED mathematical identities of a specified finite, noninteracting, lead-coupled model; CONDITIONAL device-only electrical drive and supplied projected positions; Xiang PRB original mixed orbital E/B transport remains UNRESOLVED and the manuscript audit remains PARTIAL.**

## What is now matched, and what is not

This version puts **the exact same finite-device operator used in v66** into a *different, explicitly microscopic* lead model:

\[
H(B,E)=H_B-E_bP_b(B),\qquad
J_a(B,E)=i[H(B,E),P_a(B)], \quad J_a^E=-i[P_b,P_a].
\]

The supplied Hermitian projected P_a/P_b have off-diagonal orbital dipoles and need not commute. Static magnetic B acts through the full v65 projected Peierls Hamiltonian and two lead-device hopping vectors, each rotated with the magnetic/orbital basis. The scalar potential \(-E_b P_b\) acts **only on the four-site device**, while two semi-infinite tight-binding leads remain in Fermi equilibrium. Response refers to the projected **device polarization current**, not necessarily a separately measured terminal current. In contrast, v67 measured Landauer **reservoir-voltage** transmission, which is not the same observable.

## Microscopic finite-frequency Keldysh response

Let \(\Gamma(E)=i(\Sigma^r-\Sigma^a)\), \(A(E)=G^r(E)\Gamma(E)G^a(E)\), and \(G^<(E)=if(E)A(E)\). For a device-local perturbation \(V(t)=-E_bP_b e^{-i\omega t}\) and unperturbed **fixed** reservoirs, the exact noninteracting two-energy variation is

\[
\delta G^<(\varepsilon+\omega,\varepsilon)
=-G^r_+P_bG^<_--G^<_+P_bG^a_-.
\]

Thus \(\rho_E=-i\int d\varepsilon\,\delta G^</2\pi\),
\[
\sigma_{ab}^{\rm device}(\omega,B)
=\operatorname{Tr}(J_a^0\rho_E)
+\operatorname{Tr}\bigl(-i[P_b,P_a]\rho_0\bigr),\qquad
\rho_0=\int f(\varepsilon)A(\varepsilon)\frac{d\varepsilon}{2\pi}.
\]

All \(\partial_B\) derivatives act on \(H_B\), \(P_a(B)\), \(P_b(B)\), **both retarded lead contact selfenergies**, their linewidths, and energy-shifted \(G^{r/a}_\pm\). An analytic derivative of the complete Keldysh integrand is checked independently by five-point full-field finite differences. Four response groups correspond to dB device current, dB two-energy density, derivative of explicit E-current contact, and contact times dB equilibrium density. Each term is representation dependent, the sum is not.

## A genuine two-energy reservoir Ward identity

Do not impose \([H,\rho_E]=0\), even at a site without direct electrode coupling: the electric perturbation and finite-frequency accumulation supply sources. The **matrix** continuity identity is

\[
\omega\rho_E-[H_B,\rho_E]+[P_b,\rho_0]
=\mathcal C_{\rm leads}(\omega),
\]

where \(\mathcal C_{\rm leads}\) is evaluated **independently** from two-energy Dyson and reservoir lesser selfenergies. It is not defined by subtracting the LHS. A site projector disconnected from all electrode vectors has zero projected collision; full-device trace requires nonzero lead charge injection at finite frequency. At \(\omega=.29\), the matrix Ward residual is ~9e-11, and the disconnected site residual ~1e-11, with measurable lead injection norm ~0.10547. Test remains valid under the three B-dependent representation transformations.

## Numerics and resonance-aware integration

Model: finite four-orbital v65/v66 graph, specified Peierls B=.23, q/hbar=-.84, mu=.13, T=.37, two physical single-channel semi-infinite leads with hopping 1.1, omega=.29. Units are *dimensionless chosen graph units*, not an experimental 2D conductivity. Integrate over the exact lead continuum using Gauss-Legendre quadrature split at both energy and energy+omega band edges. Narrow device resonances make low-order quadrature **qualitatively wrong**, so the reference requires 640 nodes per interval. Check 512↔640 as an independent convergence gate; 64 nodes are an explicit failing negative control. Integrated \(\int A(E)dE/2\pi\) reproduces I to ~2e-9, supporting absence of missed bound-state spectral weight **for this parameter set only**.

| Quantity, fixed lead equilibrium | Symmetric gauge | Landau and noncommuting orbital frames |
| --- | --- | --- |
| \(\sigma_{xy}^{device}(\omega)\) | \(-0.005292649616-0.001194257399i\) | same to <1e-11 |
| \(\partial_B\sigma_{xy}^{device}(\omega)\) | \(+0.007070053900-0.004146290848i\) | same to <1e-11 |
| independent 5-pt dB, step=.005 | analytic agreement <1e-10 | analytic agreement <3e-10 |

The **different-bath, same-operator** v66 scalar-gamma response is NOT numerically identical to the exact-reservoir Keldysh response. This is expected, not a source-paper discrepancy; matching an operator does not match different dissipation, boundary and field protocols.

Independent strict DC tests (omega=0) calculate finite static E in the **full** Hamiltonian \(H-E P_b\), then re-integrate the exact lead lesser density and the full contact-corrected J(E). Three-frame finite-E central differences and four-corner E/B derivatives agree with the static limit of the dynamic Keldysh result (\(\sigma_{DC}\approx0.0029963767\), \(\partial_B\sigma_{DC}\approx0.0019013505\)). The finite-model *device* equilibrium polarization-current expectation need not vanish, since explicit lead exchange cancels it in stationarity. Never confuse that device-only current with a bulk equilibrium persistent transport current.

Negative controls: deleting \(J^E\) moves the AC response by >0.04; deleting lead magnetic selfenergy derivatives or projected-position magnetic contacts breaks orbital-basis covariance and full B finite-difference identity. Positive checks: retarded analytic surface-root branch, Hermitian spectral weight positivity, total spectral completeness, omega reflection, full Dyson Ward and charge accumulation, alternative finite E/B derivatives.

## Strict remaining open problems

**This is NOT original Xiang PRB intrinsic orbital mixed E/B Hall.** We have not demonstrated (a) homogeneous field coupling to device *and* infinite crystalline reservoirs in a thermodynamically gauge-invariant way, (b) full orbital Wannier magnetic embedding and remote bands, (c) conserved macroscopic orbital current versus projected local polarization rate, (d) Ward-corrected microscopic many-body collisions/SCBA, (e) L→infinity DC/contact limit, or (f) source-SM Eq49/52 sign/factor/three-band closure. Each remains UNRESOLVED; original manuscript remains **PARTIAL**.

The v68 tests are **40 named scientific gates + one frozen-data regeneration gate**. Next v69 should prioritize a finite-lead vs semi-infinite Keldysh independent analytic/spectral oracle and a defensible same-observable boundary/gauge continuation to bulk; do not call v68 a bulk reproduction.
