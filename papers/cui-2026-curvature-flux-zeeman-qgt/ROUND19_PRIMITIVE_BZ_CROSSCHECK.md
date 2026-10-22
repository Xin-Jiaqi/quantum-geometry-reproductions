# Round 19 | Primitive reciprocal cell and Wigner-Seitz hexagon

**Scientific status:** periodic reciprocal-domain equivalence is exact; convergence of the numerical Kubo integral is UNRESOLVED. The full paper remains PARTIAL.

The two fundamental domains of a 2D reciprocal torus have equal exact normalized integrals of periodic observables:

\[
\langle\chi\rangle_{\mathrm{hex}}=
\frac{1}{(2\pi)^2}\int_{-\pi}^{\pi}dq_1\int_{-\pi}^{\pi}dq_2\,
\chi(B\mathbf q),\qquad B=
\begin{pmatrix}1&1/2\\0&\sqrt{3}/2\end{pmatrix}.
\]

Three deliberately independent sampling paths are compared using the **same
synthetic 4-band Hamiltonian, Cartesian velocity and clean cluster-interband
Kubo**:

1. Existing D6-symmetrized six-triangle Gauss rule over the first hexagonal
   Wigner-Seitz BZ, 12 order² nodes.
2. New tensor-product Gauss-Legendre rule on the reciprocal primitive
   parallelogram, order² nodes.
3. New shifted periodic midpoint grid in primitive reciprocal coordinates,
   order² nodes.

The primitive-cell output is unprojected. Projecting it onto the selected
magnetic-tensor invariant subspace would hide finite-grid symmetry errors.
Different finite quadratures need not agree; the exact BZ integral does.

## Round18 trigger and scientific warning

The previous order-24 hexagonal reference is **not an exact value**.
At order20, group 6 differs by 3.14395e-8 from order24, about 6% of the
order24 even response amplitude 5.26258e-7. For group 6m'm' the analogous
relative difference is approximately 3.9%. The symmetry projector residuals
remain around 1e-18 but do not measure numerical convergence.

## Validation / reproducibility

Eight new tests validate primitive-cell area, weights, independent sampling
constructions, analytically vanishing Fourier harmonics, reciprocal periodic
folding of hexagon points, unprojected response parity and fail-closed errors.

From the paper subdirectory run:

    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_round19.py
    python scripts/audit_round19.py --deep

**Target count: 145 tests = 137 prior + 8 new**. Only mark remote PASS
after both Python 3.11 and 3.12 workflows have finished. The audit returns
the full complex tensor, symmetry-projector residual and absolute difference
to its finite-order hexagonal comparator; it never declares automatic
convergence PASS.

This does not establish real Wannier material operators, metallic intraband
or conserving collisions, experimental/SI KME, complete Eq.(19), Table I or
a converged material response. Whole paper remains PARTIAL.


## First independent 3-method audit: actual remote CI findings

GitHub Actions run 39742758842 (Python 3.11 and 3.12) finished:
**145/145 PASS** in both jobs, including this audit. The unprojected
primitive tensor exposes significant finite-grid artifacts.

| Group | Method / size | max \|even\| | max \|even - hex12\| | raw even projection residual |
| --- | --- | ---: | ---: | ---: |
| 6 | Hexagon order12 / 1728 nodes | 4.85529e-7 | 0 | 4.60e-19 |
| 6 | Primitive Gauss20 / 400 nodes | 5.14225e-7 | 2.78444e-7 | 1.24727e-7 |
| 6 | Primitive shifted midpoint32 / 1024 nodes | 5.47477e-7 | 8.05684e-8 | 1.33997e-8 |
| 6m'm' | Hexagon order12 / 1728 nodes | 7.22893e-7 | 0 | 1.05e-18 |
| 6m'm' | Primitive Gauss20 / 400 nodes | 5.52593e-7 | 1.80528e-7 | 2.66171e-8 |
| 6m'm' | Primitive shifted midpoint32 / 1024 nodes | 6.31486e-7 | 1.18600e-7 | 1.04209e-8 |

For group 6, the primitive Gauss20 full-tensor discrepancy is **57% of
the hex12 maximum even response**, despite a superficially similar
tensor maximum magnitude. For group 6m'm', midpoint32 has
even/odd tensor deviations despite max response magnitudes that can
look similar to hex12. Comparing just tensor maximum norms is insufficient.

The primitive grid does not preserve the full point group at each finite
order. These nonzero projector residuals are **not a violation of the
Hamiltonian's symmetry**; they diagnose numerical integration. The
D6-symmetrized hexagon enforces projector zeros by construction but still
requires independent convergence testing.

The complete two-by-three complex tensors, 8 numerical orders/methods
per group, and point-count/absolute-difference evidence are frozen in
results/frozen/round19_primitive_bz_audit.json from the actual CI log.

**Scientific decision:** the method and comparison capability are VERIFIED;
BZ **numerical convergence and macroscopic response remain UNVERIFIED**.
No percentage of full-paper completion is inferred.
