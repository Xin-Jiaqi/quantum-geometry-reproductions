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
