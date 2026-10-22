# Round 17: non-Abelian Kramers response and first hexagonal BZ

**Acceptance:** limited model-method VERIFIED; the Cui paper remains **PARTIAL**.

## A. Gauge-covariant degenerate clusters

The prior single-band directed response correctly rejects exact degeneracy because a single-band position matrix element is not defined within a degenerate subspace. This round adds a separate nonabelian_response.py module rather than weakening that guard.

For a spectral cluster \(C\), the Zeeman QGT block is
\[
\mathcal T_C^{ab}=P_C\,r^a(1-P_C)\,s^b\,P_C.
\]
It is returned in a numerical basis for the cluster and transforms as
\(W_C^\dagger\mathcal T_C^{ab}W_C\); traces and eigenvalues are invariant.
The clean interband cluster Kubo response excludes intra-cluster transitions.
The independent \(N^2\) orbital Liouville solver uses the same cluster-flat
effective Hamiltonian and agrees with the spectral sum.

Acceptance tests include: exact spinful Kramers doublets, arbitrary global
unitary basis transformations, arbitrary rotations inside each degenerate
subspace, near-degenerate flattening, projector idempotency, non-Abelian QGT
trace/eigenvalue invariants, independent Kubo/Liouville equality, and the
TRIM frequency parity \(\chi_{JB}(z)=\chi_{JB}(-z)\),
\(\alpha_{ME}(z)=-\alpha_{ME}(-z)\) for the declared toy model.

**Approximation boundary:** levels grouped within the relative tolerance are
flattened to the cluster mean energy. This is exact for exact degeneracy and
an explicit cluster approximation for near degeneracy. It is NOT a conserving
metallic response, an intraband/Drude treatment, a disorder bath, or a full
non-Abelian Keldysh/Kubo-Streda theory.

## B. Triangular reciprocal lattice and the actual first hexagonal BZ

The periodic model uses reciprocal basis
\[
B=\begin{pmatrix}1&1/2\\0&\sqrt3/2\end{pmatrix},\qquad
\mathbf k=B\mathbf q,\quad \mathbf q\equiv\mathbf q+2\pi\mathbf n.
\]
The first Wigner-Seitz BZ is the regular hexagon of circumradius
\(2\pi/\sqrt3\) and area \(2\sqrt3\pi^2\). The model is periodic in the
triangular reciprocal coordinates and converts derivatives back to Cartesian
velocity components before evaluating response tensors.

The quadrature maps Gauss-Legendre nodes over six center-to-edge triangles.
The first rule was exactly \(C_6\)-symmetric but not mirror-symmetric; CI
detected a forbidden even tensor component at \(O(10^{-8})\). The corrected
rule averages the nodes and weights with one mirror image, which together
with \(C_6\) generates a fully \(D_6\)-symmetric weighted set. Its weights
sum to the exact hexagon area before normalization.

Verified in the toy model: reciprocal periodicity, integer unimodular
fractional actions for hex-compatible groups, local H/Cartesian-velocity and
full-complex-frequency Kubo covariance for \(6,6mm,6m'm',3m'\), and
integrated even/odd magnetic tensor projector checks for \(6,6m'm'\).
The order-3 integrals are **symmetry checks only**, not converged physical
predictions. Incompatible square-lattice groups such as \(4mm\) are rejected.

## CI and reproducibility

- 19 new tests; 137 total.
- Python 3.11 and 3.12 each report **137 tests, OK** after the mirror-quadrature fix.
- Reproduce from the paper directory:
      python -m pip install -e .
      python -m unittest discover -s tests -v
      python scripts/audit_round17.py

## Still unresolved

No real Wannier/material response, physical spin/position operator interface,
conserving collision or contact terms, non-Abelian intraband metallic KME,
converged BZ observable, full Table I adjudication, or SI conversion. The
method validates bounded microscopic model identities only.
