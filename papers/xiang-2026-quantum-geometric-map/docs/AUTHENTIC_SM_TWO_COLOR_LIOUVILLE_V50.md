# v50 | Authentic SM two-frequency Liouville: Eq.(4) to Eq.(9) to Eq.(20) versus Eq.(31)

**Primary source:** the user's genuine 8-page PRB 113 L201406 and 16-page Supplemental Material (SM), directly inspected and SHA256-pinned. Original PDFs not copied into GitHub. PRB SHA256: b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730. SM SHA256: 9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4. **Overall original-paper acceptance PARTIAL.**

## 1. Exact original time-dependent two-colour perturbation protocol

SM pp.3-4 Eqs.(4)–(10) use separate omega_a,omega_b, eta_a,eta_b with two field orderings. The magnetic perturbation uses **Bbar = -g mu_B B**. In a *fixed physical-spin basis* define

\[
L(z)X=zX+i[H_0,X],\qquad z_a=\eta_a-i\omega_a,\quad z_b=\eta_b-i\omega_b.
\]

Then the original SM density hierarchy, with one Fourier amplitude per field, is

\[
\rho_E=L(z_a)^{-1}\partial_{k_b}f(H_0),\qquad
\rho_{\bar B}=L(z_b)^{-1}\{-i[\sigma_x,f(H_0)]\},
\]
\[
\boxed{\rho_{E\bar B}=L(z_a+z_b)^{-1}
\left(\partial_{k_b}\rho_{\bar B}-i[\sigma_x,\rho_E]\right)}.
\]

The Fourier components are divided by 2 in the source's explicitly real E(t), Bbar(t); this implementation outputs a response **per the product of individual complex Fourier-component amplitudes**. It does not conflate that amplitude convention with published physical Fig2 units.

The solver checks first-order **original SM Eq.(6) and (7)** against their eigenbasis matrix fractions at nonzero omega and eta. A separate matrix-identity residual validates inversion of L. The entire f(H0) density matrix and rho_B are differentiated with respect to k; no manual Berry-connection phase or v46 analytic source shortcut is used.

**SM p.5 Eq.(15)** introduces a phenomenological tau **only to regularize a selected divergent 1/eta_a electric adiabatic pole**. It does not add a finite uniform eta to every coherence channel. v47/v48 frozen-H0 constant-eta RTA is a separate method, and the source-selected pole is obtained here without silently imposing that different collision model.

## 2. Independent source Eq.(9) solution converges to the sign of SM Eq.(20)

For main Eq.(10) physical spin helicity, H=beta*kx+kx*sigma_y-ky*sigma_x, beta=1/8, mu=.3,T=.03, radial Gauss100 by angular128, cutoff2.5, e=hbar=vF=gmuB=tau=1:

| eta_a | eta_b | eta_a x full two-field current divided by physical B |
| ---: | ---: | ---: |
| .010 | .0005 | **-0.00250733894** |
| .002 | .002 | -0.00250654116 |
| .002 | .0002 | -0.00250624973 |
| .0002 | .002 | -0.00250624965 |
| .0005 | .0005 | **-0.00250621056** |
| Direct original SM Eq.(20) | | **-0.00250618825** |
| **Printed SM Eq.(31) on same helicity** | | **+0.00250618825** |

The **unmodified original finite-two-frequency Eq.(4)/(9)** independently agrees with the sign of source Eq.(20), *not* the printed Eq.(31) sign. This does not establish the measured Hall-voltage sign, because an explicit physical bath protocol and contact contributions remain unproved.

## 3. Field ordering and AC/DC limits

At eta_a=.002,eta_b=.001, the eta_a-weighted **B then E** term is +0.004994807, and **E then B** is -0.007501182; their correct total is **-0.002506375**. Changing the eta_a/eta_b ratio strongly changes the separate contributions; it leaves their combined small-regulator dc pole close to the original Eq.(20) across tested regulator paths. Neither field ordering should be used alone to accept source Eq.(31).

With omega_a=.04, omega_b=0 at the same eta, eta_a times the EB current is approximately **-6.46e-6 - 1.25e-4 i**, *not* the -0.002506 static pole. Thus a fixed nonzero electric frequency removes the particular 1/eta_a DC singularity. With omega_b=.05 the complex response also differs from DC; (omega_a,omega_b) -> (-omega_a,-omega_b) gives the complex conjugate, independently tested.

## 4. Scientific status and remaining gates

**PASS scientific method subset:** original finite frequency Eq.(4)/(6)/(7)/(9) in physical-spin matrices; correct magnetic-field sign; both electric/magnetic orderings; original selected tau pole Eq.(15)/(20); two tensor directions and both main/SM helicities; independent eta regulator paths, finite-frequency conjugation, k finite-difference convergence, source sign and input validation. 12 new regression tests.

**OPEN:** printed SM Eq.(20) -> Eq.(31) relative sign; main Eq.(10) versus SM Eq.(62) opposite spin helicity and transformed physical sigma_x; complete steady-state bath/collision interpretation; original Fig.2(b) absolute SI ordinate roughly 2pi away from literal expression; orbital higher-band and three-band processes. **CI PASS does not prove full original-paper physics; acceptance remains PARTIAL.**

Reproduce inside the isolated PRB module:

    python -m unittest tests.test_sm_two_color_liouville_v50 -v
    python scripts/reproduce_sm_two_color_v50.py --eta-a .002 --eta-b .001 --nr 100 --nphi 128

For the archived independent regulator/AC scan, see results/frozen/v50_two_color_regulator_frequency_scan.json. Optional original PDF paths --sm-pdf and --main-pdf verify exact source bytes by SHA256. Raw copyrighted PDF bytes are not uploaded.

**v51 priority:** independent source index/velocity and position-operator sign audit for Eq.(31) and the complete spin plus B transformation, followed by a signed, **unfitted** absolute Fig.2(b) SI analysis. Whole paper remains PARTIAL.
