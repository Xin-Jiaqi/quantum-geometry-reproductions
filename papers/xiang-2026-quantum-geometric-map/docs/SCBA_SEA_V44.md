# v44 | SCBA energy-dependent Hall sea II and provisional I(a)+II sum

**Scientific state: PARTIAL.** This is an independently constructed **linear charge-Hall QWZ method**, not Xiang et al. PRB 113 L201406's nonlinear mixed EB PHE, nor a proof of a conserving disorder-averaged Středa conductivity. The package retains an exact finite-torus v42 comparison and explicitly reports disagreement.

## One microscopic disorder convention

Use the QWZ Hamiltonian with m=-1 and physical lattice velocities v_a=∂H/∂k_a, chemical potential mu=1.1, bare Green-function regulator eta=0.12, and Gaussian scalar onsite disorder **variance W** in both SCBA and exact torus. No independent fitting of scattering times. Response units e²/h.

For every real integration energy E, self-consistently solve

\[
G^R(E,k)=[E+i\eta-H(k)-\Sigma^R(E)]^{-1},\quad
\Sigma^R(E)=W\langle G^R(E,k)\rangle_{\mathrm{BZ}},\quad G^A=(G^R)^\dagger.
\]

The matching RR **density** vertex obeys

\[
\Gamma_0^{RR}=I+W\langle G^R\Gamma_0^{RR}G^R\rangle,\quad
\partial_E\Sigma^R=I-\Gamma_0^{RR},\quad
\partial_EG^R=-G^R\Gamma_0^{RR}G^R.
\]

This is verified independently using freshly converged Sigma(E+h), Sigma(E-h), G(E+h), G(E-h). A naive -G² misses the energy dependence of the self-energy. For a momentum-independent SCBA Sigma on a complete periodic BZ the **RR spatial** current-ladder correction vanishes by periodic momentum Ward, which motivates bare spatial velocities in the sea below. This is a limited model identity and **not** full electromagnetic Ward/contact closure.

## Explicit response pieces

Reuse v43's RA current ladder on **one** side of the antisymmetric Středa I(a) surface:

\[
\sigma_H^{I(a)}/(e^2/h) = \tfrac12\mathrm{Re}\langle\mathrm{tr}
[v_xG^R\Gamma_y^{RA}G^A-v_yG^R\Gamma_x^{RA}G^A]\rangle_k.
\]

At each energy quadrature node, solve Sigma and Gamma0 and compute the sea candidate

\[
\sigma_H^{II,\mathrm{candidate}}/(e^2/h)
=\tfrac12\int_{-\infty}^{\mu}dE\,\mathrm{Re}\langle\mathrm{tr}
[v_xG^Rv_y\partial_EG^R-v_x\partial_EG^Rv_yG^R-(R\to A)]\rangle_k.
\]

The API intentionally returns **candidate_total**, not physical_hall_conductivity. The conventional surface/sea split contains overlap and is not a unique physical mechanism partition. Both parts are computed under the **same** W and eta.

## Independent checks and numerical convergence

- Clean W=0: directly integrated sea plus surface versus independently diagonalized two-band Bloch pole primitive; mesh 35², 420 quadrature nodes from -12 to mu, candidate total -0.2928979048, clean analytic -0.2928990006 e²/h, residual 1.10e-6 dominated by finite cutoff.
- Finite W=0.03 and 0.12: RR density Ward/BSE residual and **independently** finite-differenced energy derivatives agree at E=-2.1, -0.75, 0.82. The energy-integrand derivative is independently checked without the ladder formula.
- A coarse global energy grid (80–160 nodes) produces unacceptable errors in some clean metallic cases, even when SCBA fixed-point residuals are excellent. Frozen reference uses **420 Gauss–Legendre nodes**, lower energy -12, with an independent 360-node/-9 convergence check. Smaller eta requires stronger integration refinement.
- Independent v42 exact reference: full eigenstate Středa I(a)+II in frozen real-space disorder. Twist-average within each sample; sample SEM is computed **only across independent disorder seeds**, not twists.

## Frozen results

m=-1, mu=1.1, eta=0.12. SCBA mesh 35×35; exact L=4, 24 independent disorder seeds 1947–1970, 7×7 twists each; e²/h.

| W | Approx. I(a) | Approx. II | **SCBA candidate sum** | Exact I(a) | Exact II | **Exact sum ± seed SEM** | Candidate minus exact |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.03 | -1.810174 | +1.434051 | **-0.376124** | -1.854201 | +1.464349 | **-0.389852 ± 0.021380** | +0.013729 |
| 0.12 | -1.430324 | +0.974264 | **-0.456059** | -1.500275 | +1.016586 | **-0.483690 ± 0.027788** | +0.027630 |

**Interpretation:** the combined approximation is close to the finite-sample exact result in these defined cases. A small difference relative to the sampling SEM does not bound systematic errors from finite size, crossed disorder diagrams, or missing electromagnetic Ward and contact terms. It cannot establish a full conserving SCBA theory or the original paper's nonlinear EB PHE.

## Reproduce

~~~bash
cd papers/xiang-2026-quantum-geometric-map
python -m pip install -e .
python -m unittest tests.test_scba_sea_v44 -v
python scripts/reproduce_scba_sea_v44.py --W .03 --mesh 35 --L 4 --seeds 24 --ntheta 7 --lower -12 --quadrature 420
python scripts/reproduce_scba_sea_v44.py --W .12 --mesh 35 --L 4 --seeds 24 --ntheta 7 --lower -12 --quadrature 420
~~~

## Remaining scientific gate

v45 needs an explicitly conserving response functional, all required transverse electromagnetic Ward/contact identities, and systematic size/disorder/eta/DC limit checks. Only after a genuine transport-level derivation should a claim of complete averaged Středa Hall be evaluated. The original PRB steady-state mixed EB response, collision/bath protocol, multiband orbital terms, and Fig.2(b) normalization discrepancy remain **OPEN**.
