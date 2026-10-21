# v42: Exact finite-sample Smrčka–Středa decomposition against v41 Hall

**Paper being independently studied:** Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB **113**, L201406 (2026). **Overall scientific status: PARTIAL.** This v42 is an independent *linear* charge-Hall test on a QWZ model; it is **not a formula claimed from the PRB**, an implementation of the original second-order mixed-field \(EB\) response, or a completed conserving infinite-system SCBA–Středa calculation.

## 1. One microscopic Hamiltonian and observable

Reuse the exact **same** `qgmap.exact_hall_v41.torus` for both Hall approaches:

\[
H_{\theta}=F\,\mathrm{diag}_{\bf k}\big[\sin(k_x+\theta_x/L)\sigma_x+\sin(k_y+\theta_y/L)\sigma_y+[m+\cos(k_x+\theta_x/L)+\cos(k_y+\theta_y/L)]\sigma_z\big]F^\dagger+V.
\]

The site-diagonal scalar disorder \(V\) has variance \(W\) and is **frozen** when the boundary twist changes. Physical single-particle charge-current vertices in this protocol are \(J_a=\partial H/\partial A_a\), with \(A_a=\theta_a/L\). Use fixed \(T=0\), \(e=\hbar=a=1\), chemical potential \(\mu\), and a common positive scalar Green regulator \(\eta\):

\[
G^{R,A}(\epsilon)=(\epsilon-H\pm i\eta)^{-1},\qquad
\partial_\epsilon G^{R,A}=-(G^{R,A})^2.
\]

The antisymmetric Hall channel is \(\sigma_H=(\sigma_{xy}-\sigma_{yx})/2\) in \(e^2/h\), with the sign convention fixed in v38/v41. For this QWZ model with an independent scalar onsite \(V\), the mixed **linear Hall current contact** \(\partial^2H/\partial A_x\partial A_y\) is exactly zero, verified by v41 finite differences. This statement must **not** be exported to another Hamiltonian or to the nonlinear \(EB\) response.

## 2. Explicit standard Smrčka–Středa I/II (not a unique mechanism decomposition)

Let \(N=L^2\), and \(\mathrm{Tr}\) denote the complete \(2N\)-dimensional one-particle trace. Define at \(\epsilon=\mu\)

\[
\begin{aligned}
\frac{\sigma^{I(a)}_{xy}}{e^2/h}&=\frac1N\mathrm{Re\,Tr}[J_xG^RJ_yG^A],\\
\frac{\sigma^{I(b)}_{xy}}{e^2/h}&=-\frac1{2N}\mathrm{Re\,Tr}[J_xG^RJ_yG^R+J_xG^AJ_yG^A],\\
\frac{\sigma^{II}_{xy}}{e^2/h}&=\frac1{2N}\int_{-\infty}^{\mu}d\epsilon\,\mathrm{Re\,Tr}\big[J_xG^RJ_y\partial_\epsilon G^R-J_x\partial_\epsilon G^RJ_yG^R-(R\to A)\big].
\end{aligned}
\]

We retain **only the antisymmetrized Hall part** of each term. For the same Hermitian \(H,J_x,J_y\), \(\sigma^{I(b)}_H=0\) by cyclicity; \(I(a)_H\) is **not** zero in general. This differs from v40's isolated antisymmetric **RA correlation kernel**, whose overall physical normalization and complementary sea term were not fixed there. The two cannot be numerically equated by renaming them.

References for this **standard response decomposition**, distinct from the Xiang paper: Sinitsyn et al., [PRB **75**, 045315 (2007)](https://doi.org/10.1103/PhysRevB.75.045315). The **surface/sea assignment is not unique**; conventional Smrčka–Středa has an overlap under integration by parts, as analyzed by Bonbien and Manchon, [PRB **102**, 085113 (2020)](https://doi.org/10.1103/PhysRevB.102.085113). Hence \(I\) and \(II\) are not independently measured "intrinsic" and "extrinsic" conductivities. The **sum**, after the declared limits, is the quantity compared to v41.

## 3. Pairwise energy integral and independent code paths

Diagonalize the *entire finite disordered torus*, obtaining \(E_n\) and \(P_{nm}=(J_x)_{nm}(J_y)_{mn}\). For each \(n<m\) write \(\Delta=E_m-E_n>0\), \(A_{nm}=\mathrm{Im}P_{nm}\), \(z_n=\mu-E_n+i\eta\) and \(z_m=\mu-E_m+i\eta\). The v38 two-pole primitive generalizes from an isolated \(2\times2\) Bloch block to **all eigenstate pairs**:

\[
\begin{aligned}
F_{nm}&=-\frac{1}{\Delta}\left(\frac1{z_n}+\frac1{z_m}\right)+\frac{2}{\Delta^2}(\log z_n-\log z_m),\\
\frac{\sigma_H^{II}}{e^2/h}&=-\frac2N\sum_{n<m} A_{nm}\,\mathrm{Im}F_{nm},\\
\frac{\sigma_H^{I(a)}}{e^2/h}&=\frac2N\sum_{n<m}\frac{A_{nm}\eta\Delta}{[(\mu-E_n)^2+\eta^2][(\mu-E_m)^2+\eta^2]}.
\end{aligned}
\]

Independent checks deliberately avoid the above pairwise formula: `direct_surface_trace` uses a **dense inverse-matrix trace** at \(\mu\); `direct_sea_quadrature` integrates the full \(G^{R,A}\) matrix expression over **real energy** with \(\partial_\epsilon G=-G^2\), *without* eigenvector primitives. The lower energy cutoff is explicit and checked. In the clean limit, the same Bloch mesh reproduces the existing v38 sea integral. Phase/sign and basis checks include switching \(x\leftrightarrow y\) and an arbitrary constant unitary Hilbert-basis transformation.

**Critical regulator distinction:** \(\eta\) here is a chosen one-body Green regulator. The independent v41 \(\omega=0,\delta=0\) Lehmann result is recovered as \(\eta\to0^+\) at **fixed finite sample and a chemical potential separated from every eigenlevel**. Equality is NOT asserted for a generic finite \(\eta\), and the order of the thermodynamic/disorder/frequency limits has not yet been settled. Neither \(\eta\) here nor the v39 external Green regulator can be identified with a collision rate without deriving a conserving scattering approximation.

## 4. Reproduced benchmarks (units \(e^2/h\))

Use \(L=4,m=-1,W=0.12,\mathrm{seed}=1947,\boldsymbol\theta=(0.11,0.33)\), *the same one frozen onsite realization*.

| \(\mu\) | \(\eta\) | Hall I(a), surface | Hall II, sea | **Total** | independent v41 Lehmann | |total−v41| |
|--:|--:|--:|--:|--:|--:|--:|
| 0 | .003 | −.005593286 | −1.269821006 | **−1.275414292** | −1.275414375 | 8.33e−8 |
| 1.1 | .16 | −1.101567325 | +.901654789 | **−.199912536** | −.307005888 | .10709 |
| 1.1 | .05 | −1.188019579 | +.913215253 | **−.274804326** | −.307005888 | .03220 |
| 1.1 | .01 | −.376276279 | +.069713870 | **−.306562410** | −.307005888 | 4.43e−4 |
| 1.1 | .003 | −.115816488 | −.191177158 | **−.306993647** | −.307005888 | 1.22e−5 |

For the \(\mu=1.1\) bulk-metallic parameter, the **finite torus has discrete levels**: the nearest level is about 0.04970 model energy units from \(\mu\). The large and even sign-changing individual terms at finite \(\eta\) must not be interpreted as a material's measured Fermi-surface/sea partition. A six-seed **finite** disorder sample mean is additionally recorded, not extrapolated to the thermodynamic limit and **not** used as an SCBA acceptance.

## 5. Reproduction and scientific acceptance

~~~bash
cd papers/xiang-2026-quantum-geometric-map
python -m pip install -e .
python -m unittest tests.test_exact_streda_v42 -v
python scripts/reproduce_exact_streda_v42.py --L 4 --mass -1 --W .12 --mu 1.1
python scripts/reproduce_exact_streda_v42.py --L 4 --mass -1 --W .12 --mu 1.1 --etas .005 --ensemble-seeds 1947 1948 1949 1950 1951 1952
~~~

Ten named v42 tests cover exact clean/disordered Lehmann recovery, metallic cancellation, independent direct matrix trace, independent quadrature, equivalence to v38 Bloch sea, unitary basis / Hall sign, a small fixed seed ensemble, convergence, and fail-closed inputs.

**Physics gate PASS:** finite-size one-body *linear* Hall standard Středa I+II reproduces full finite-system Kubo–Lehmann as \(\eta\to0\) on these gapped finite samples and bulk-metal parameter (which is also discrete at finite \(L\)). **Still OPEN:** conserving *disorder-averaged* SCBA/Kubo–Bastin I+II with energy-dependent self-energy and matched vertex classes, disorder and size scaling, actual metallic thermodynamic DC, external-field collision/bath preparation, original PRB nonlinear \(EB\) PHE and multiband orbital corrections. **Do not upgrade whole-paper status from PARTIAL.**
