# v40: Matched SCBA current ladder and antisymmetric-Hall zero controls

**Scientific state: PARTIAL.** This calculation is an independent QWZ two-band method test, **not** a calculation of the full Xiang et al. PRB 113 L201406 nonlinear DC response. All quantities are dimensionless Brillouin-zone means; an isolated RA bubble cannot be converted into a sheet Hall conductivity.

## Definitions

Use the v39 *same energy-dependent scalar disorder SCBA* in the QWZ lattice model,

\[
G^{R}(E,k)=[E+i\eta-H(k)-\Sigma^R(E)]^{-1}, \quad
\Sigma^R(E)=W\langle G^R(E,k)\rangle_k,
\]

and \(G^A=(G^R)^\dagger\). For \(s,t\in\{R,A\}\) each momentum-dependent current vertex is

\[
\Gamma_a^{st}(k)=v_a(k)+C_a^{st},\qquad
C_a^{st}=W\langle G^s(k)[v_a(k)+C_a^{st}]G^t(k)\rangle_k.
\]

The k-independent \(C_a^{st}\) is solved using an independent four-component Pauli basis, in each of the RA, AR, RR and AA channels, with **the same variance \(W\)** as the self-energy. No separate relaxation time is fitted.

The fixed-energy channel kernel is the BZ mean

\[
B_{ab}^{st}=\langle\mathrm{Tr}[v_aG^s\Gamma_b^{st}G^t]\rangle_k.
\]

A symmetric spectral diagnostic is assembled as

\[
S_{ab}=\frac14(B_{ab}^{RA}+B_{ab}^{AR}-B_{ab}^{RR}-B_{ab}^{AA}).
\]

This is the specific four-channel **spectral** kernel already used in v18. It is **not** the complete Kubo–Bastin or Středa conductivity. In particular, its antisymmetric part vanishes (independently verified for this noncrossing model), even though \(B_{xy}^{RA}-B_{yx}^{RA}\) is nonzero.

## Independent Ward checks

**RA density vertex.** The *same* SCBA \(W\) gives, at fixed regulator \(\eta>0\),

\[
\Gamma_0^{RA}=I+W\langle G^R\Gamma_0^{RA}G^A\rangle_k
=I-\frac{\Sigma^R-\Sigma^A}{2i\eta}.
\]

This is a **resolvent identity**, not a demonstration that a physical finite-scattering-rate charge current obeys continuity.

**RR current vertex.** In a translationally invariant periodic BZ with momentum-independent self-energy,

\[
\partial_{k_a}G^R=G^Rv_aG^R,\quad
\langle G^Rv_aG^R\rangle_k=\langle\partial_{k_a}G^R\rangle_k=0.
\]

Thus the continuum-BZ RR current ladder correction must vanish. Its small **discrete-grid** remnant converges to zero on dense meshes and must not be confused with a finite physical current renormalization.

## Frozen numerical benchmarks

QWZ \(m=-1\), \(\mu=1.1\), \(W=0.12\), external \(\eta=0.12\), \(N=21\):

| Calculation | Dimensionless BZ kernel |
| --- | ---: |
| Bare RA antisymmetry \((B_{xy}^{RA}-B_{yx}^{RA})/2\) | −1.27886213285 |
| Current-ladder dressed RA antisymmetry | −1.43106823404 |
| **RA vertex effect**, dressed minus bare | **−0.15220610119** |
| Four-channel spectral antisymmetry, absolute | \(2.95\times10^{-17}\) |
| RA density Ward residual | \(1.28\times10^{-11}\) |

RR current-ladder correction on grids \(N=21,31,41,61\) decreases from approximately \(3.82\times10^{-4}\) to \(2.04\times10^{-9}\). For \(W=0\), every ladder correction vanishes identically.

**Independent real-space torus cross-check.** An explicit \(L=4\) Fourier-transformed QWZ Hamiltonian and full single-particle diagonalization give the bare RA antisymmetric kernel **−1.93552124267** at \(W=0\), identical to the unshifted four-by-four clean Bloch mesh. With \(W=0.12\), six fixed disorder realizations give a sample mean near **−2.01220**, with standard deviation near **0.45350**. These samples are **not an infinite-system SCBA benchmark**: neither size nor ensemble averaging is converged, and a disorder-dressed current vertex is absent from the exact finite-sample *bare* observable.

## Scientific conclusions and remaining work

1. The matching scalar-disorder RA current ladder can materially change an *isolated* antisymmetric RA correlator in a metallic QWZ model. **That is not yet a physical Hall current.**
2. The symmetric four-channel spectral sector cancels antisymmetric Hall contributions in its own response definition; consequently it cannot be substituted for a complete Hall Kubo–Středa/Bastin calculation.
3. v39's RR **density** vertex, this round's RA **current** vertex and v38's clean Fermi-sea integral are different response objects. Simply summing them is not a controlled derivation of a physical DC conductivity.
4. A complete finite-energy disorder-dressed Hall calculation needs consistent Fermi-surface plus Fermi-sea terms, all required vertices and contact terms, and separately controlled \(q\to0\), \(\omega\to0\) and bath broadening limits. Crossed diagrams may be needed beyond noncrossing SCBA. Original paper's second-order mixed \(EB\) PHE is additionally nonlinear and not addressed here.

## Run

```bash
python -m pip install -e .
python -m unittest discover -s tests -p 'test_current_vertex_v40.py' -v
python scripts/reproduce_current_vertex_v40.py --mesh 21 --energy 1.1 --W .12 --eta .12
python scripts/reproduce_advanced.py --case scba-current-vertex-v40
```

The module validates Bethe–Salpeter residuals, density Ward, RR momentum Ward with mesh convergence, SU(2) basis invariance, clean finite-torus equivalence, and explicit zero-Hall control. A passing CI only certifies these **named model methods**.
