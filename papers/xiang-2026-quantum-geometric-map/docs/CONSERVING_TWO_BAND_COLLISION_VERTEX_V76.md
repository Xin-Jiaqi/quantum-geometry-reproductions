# v76 | Two-band elastic collision, coherent GKSL, and an explicit magnetic collision vertex

**23 October 2026. Original article:** Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113 L201406 (2026). **Whole-paper science status: PARTIAL.** Sources: this repository's previously verified v49/v50/v53 SM audits, and the frozen v75 handoff. v76 provides independently defined model results, **not source-equivalent spin PHE or orbital MNHE**.

## What genuinely closes

### I. Both bands and multiple true isoenergy shells for the article's spin Bx Dirac geometry

Take the physically specified continuum tilted Dirac Hamiltonian (model sign convention)

\[
H(k,B_x)=\beta k_x I+k_x\sigma_y-(k_y+B_x)\sigma_x,\quad |\beta|<1.
\]

For positive energy magnitude \(E\) and band \(b=\pm1\), **the exact elastic shell**, written in terms of the spinor angle \(\theta\), is

\[
r_b(E,\theta)=\frac{E}{1+b\beta\cos\theta},\quad
k_x=r_b\cos\theta,\quad k_y=r_b\sin\theta-B_x,\quad
\varepsilon_b=bE.
\]

The **band-specific DOS weights** are proportional to \((1+b\beta\cos\theta)^{-2}\). Use a scalar-impurity Born-style kernel \(W_{i\to j}^{(b,E)}=\gamma w_{j,b}[1+\cos(\theta_i-\theta_j)]\) for \(i\ne j\); same-band and same-energy scattering only. Both bands have detailed balance, conservation of each energy-shell occupation, positive relaxation spectrum, and direct-versus-independent weighted spectral current oracles. The transport current is integrated over **two explicit energy magnitudes** and the Fermi-derivative factors for both bands. Results are dimensionless and omit physical continuum integration normalization; this is an **energy quadrature demonstration**, not a completed transport coefficient in SI units.

With E=(0.3,0.6), beta=.125, gamma=.07, mu=.19 and temperature=.16, the two-band modeled current is **10.962386926061** in the declared normalization, with +band 9.451833734794 and -band 1.510553191267. Rigid Zeeman shift Bx changes ky of the shells but not their intrinsic spinor angles, velocities or short-range impurity transfer q, so this particular **population-only** DC and AC response has \(\partial_{B_x}J=0\). That is **not** a no-go proof for the article's *interband-coherent* spin PHE. Negative control: holding a fixed laboratory-radius ring at Bx=.08 yields an energy spread ~0.17571 and violates elastic scattering; it must be rejected.

### II. A charge-conserving completely positive two-band Liouville dissipator

To represent **interband coherence at identical momentum**, independently restrict to the untilted beta=0 geometry: at each theta, the + and - bands have the same physical radius. Use the gapped Dirac family

\[
H_0(\theta)=E(\cos\theta\,\sigma_y-\sin\theta\,\sigma_x)+m\sigma_z.
\]

Band populations \(p_{\theta,b}\) are coupled by positive elastic spinor-overlap jumps \(L_{j\leftarrow i,b}=\sqrt{W_{i\to j,b}}\,|j,b\rangle\langle i,b|\), with angular overlap \(2|\langle u_{i,b}|u_{j,b}\rangle|^2=1+(E^2\cos\Delta\theta+m^2)/(E^2+m^2)\). Each local interband coherence suffers the **Gorini–Kossakowski–Sudarshan–Lindblad (GKSL)** outgoing-rate dephasing \((\Gamma_{+,i}+\Gamma_{-,i})/2\), plus an independently selectable nonnegative pure-dephasing rate \(\gamma_\phi\). This block-diagonal-in-momentum GKSL restriction is charge conserving, maintains Hermiticity and positive states under exact time propagation, and has \(\mathcal D_0[f(H_0)]=0\).

For controlled *local interband vertex* \(r_y^{+-}=v_y^{+-}/(i(\varepsilon_+-\varepsilon_-))\) and spin \(S_x=\sigma_x\), compute **both** mixed Liouville orders with distinct electric and magnetic frequencies. Two independent algorithms, angular-population spectral inverse and full dense orbital superoperator, agree; a nontrivial angle-dependent U(1) band gauge does not change observables. The angularly **integrated** candidate mixed current vanishes by this model's symmetry, while its angular first harmonic is **-0.056942446+0.800529766i** (local diagnostic, not a bulk source PHE). **Crucially, the local vertex excludes the full length-gauge diagonal electric drift**, so no author's selected 1/eta electric pole is claimed here.

### III. Strong negative result: a measured transport lifetime does not determine coherent mixed response

Adding \(\gamma_\phi=0.14\) leaves the **population collision matrix and its every transport eigenvalue exactly unchanged**, but changes local mixed response and the uniform spin susceptibility. The first-harmonic mixed difference is nonzero (~0.213 in absolute magnitude); the spin susceptibility changes by >0.05. Thus the full two-band response cannot be reconstructed from the **single scalar transport tau** in the absence of coherence/dephasing information. This is a constructive GKSL counterexample, not just a heuristic argument.

### IV. Collision-field Ward: an actual B-dependent bath, not an after-the-fact numerical cancellation

For an equilibrium-prepared magnetic field, \(F_B=f_T(H_B-\mu)\) and a microscopically admissible field-dependent collision \(\mathcal D_B\), stationarity requires \(\mathcal D_B[F_B]=0\). Differentiating yields the **necessary collision contact identity**

\[
\boxed{(\partial_B\mathcal D_B)_0 F_0
=-\mathcal D_0\big[\partial_B F_B|_0\big].}
\]

For the original-orientation **Bx** Zeeman perturbation at fixed laboratory shell, freezing \(\mathcal D_0\) fails this equilibrium Ward: the local test residual is **0.3467493499** for specified model parameters; the required collision-contact trace is ~0. This is a *necessary-condition audit*, **not** a physically constructed Bx collision operator; its actual energy shell must shift as in Part I.

To show sufficiency of a physically constructed field-dependent bath rather than inserting a cancellation term, separately construct an explicit fully elastic scalar-impurity GKSL family for **Bz mass-like Zeeman** \(H_B=E(\cos\theta\sigma_y-\sin\theta\sigma_x)+(m+B_z)\sigma_z\). **This is not the article's Bx geometry.** All angles at each Bz have the same two band energies, and the jump operators and spinor-overlap rates are rebuilt at each Bz. Accordingly \(\mathcal D_{B_z}[F_{B_z}]=0\) by direct evaluation, and differentiating the **full** Liouvillian includes magnetic current commutator and \(\partial_{B_z}\mathcal D\). A fixed physical orbital probe \(V_E=\sigma_x\), measured physical \(J=\sigma_y\) (a specified finite-dimensional optical test, **not homogeneous physical electric length gauge**) yields:

| Bz derivative protocol | Complex result | Discrepancy vs independently re-diagonalized 5-point |
| --- | --- | --- |
| Full physical field derivative | \(0.356798105221-0.254547840247i\) | **3.71e-11** |
| Independent full-field 5-point | \(0.356798105242-0.254547840278i\) | reference |
| **Incorrect:** freeze collision vertex | \(0.277554354631-0.262946504900i\) | **0.0796876** |

The Bz equilibrium-collision Ward defect is **1.54e-10**, and the undisplaced B=0 state is stationary within ~1e-16. This proves the role of the collision magnetic derivative in a **fully specified conserving microscopic toy bath**. The selected \(B_z\) example cannot be passed off as the original \(B_x\) spin PHE.

## Correct scientific classification

- **VERIFIED:** true two-band/energy-shell elastic population transport; full angular DOS and Bx shift; detailed balance and U(1); charge-conserving completely positive two-band GKSL; independent dense/spectral two-color; coherent lifetime counterexample; necessary Bx equilibrium contact; **explicit** source-distinct Bz collision-Ward and independent derivative.
- **CONDITIONAL:** positive scalar short-range Born-style impurity model; specified finite angular/energy grids; coherence-only local interband electric vertex; prescribed dephasing; Bz optical probe chosen to isolate the collision derivative.
- **UNRESOLVED:** original author's full homogeneous \(B_x\) two-color covariant length-gauge electric derivative **and field-dependent conserving collision** on moving exact energy shells; original SM Eq20/31 irreducible relative sign and Eq10/62 helicity; exact steady-state tau correspondence; orbital MNHE full Wannier magnetic positional shift, Eq49/52 coefficient and omitted three-band terms; published Fig2 2pi ordinate normalization and SI values. Whole paper **PARTIAL**.

## Reproducibility

`src/qgmap/conserving_two_band_collision_v76.py` (public); `tests/test_conserving_two_band_collision_v76.py` (**52 named independent gates**); `tests/test_conserving_two_band_collision_v76_frozen.py` (**one independent frozen rerun**); `scripts/reproduce_v76.py`; `results/frozen/v76_conserving_two_band_collision.json`. Local **53/53 PASS**, with full GitHub CI and independent merged main as separate release requirements.
