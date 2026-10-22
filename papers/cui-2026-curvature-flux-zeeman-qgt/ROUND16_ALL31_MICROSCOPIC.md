# Round 16: all 31 MPG generator-to-microscopic bridge

Source: Cui et al. (2026), arXiv:2604.09725v3, Eqs.(16)–(19), Table I;
pinned 31-group data in data/toolkit_2d_generators.json.

**Acceptance: local clean interband model-method VERIFIED; full paper PARTIAL.**

## Essential corrected physical convention

The pin-set 2D magnetic generator's 2x2 **R is a real-space operation**.
Under time reversal, crystal momentum changes sign independently:
  Q=(-1)^anti R, k'=Qk.
A physical real-spin axial 3-vector transforms as
  L=(-1)^anti det(R) diag(R,1).
A=U K^anti is the projective spinor implementation. For the declared
nondegenerate clean *interband-only* Kubo protocol, at any complex z,
  chi_JB(Qk,z)=Q chi_JB(k,(-1)^anti z) L^T
  alpha_ME(Qk,z)=(-1)^anti L alpha_ME(k,(-1)^anti z) Q^T.

Not using Q (e.g. incorrectly setting time reversal k'=k because R=I)
breaks the gray-group physics. The finite-eta antiunitary branch uses
-z with *negative* imaginary part, NOT -omega+i eta.

## Independent numerical acceptance

1. All 31 pinned magnetic point groups: explicit spin-1/2 projective
   generator closure, physical-spin axial/TR parity, 4-band orbitally
   mixed Hamiltonian and analytic velocity covariance. At the generic
   k=(0.41,-0.32), z=0.051+0.011i, mu=0.17, T=0.17 all 31 pointwise
   J/B and spin-M/E spectral Kubo constraints hold. Global maximum
   operator or response residual **9.1552e-16** (natural units).
2. Selected 7 representative groups: completely independent 16x16
   orbital Liouville solution vs the same microscopic group constraints.
3. All 31 x 2 = **62 independent** even/odd in-plane J/B rank checks
   match the previous 31-group static real Jahn-character tensor projection.
   The exact bridge is Q and L, not an assumption that antiunitary operations
   leave k unchanged. This is a statement about allowable tensor spaces.
4. Symmetry-constrained clean square-BZ 8x8 midpoint examples for
   11',2',m',4',4mm: full Kubo even/odd response matches group projectors.
   For 11', the odd term cancels to roundoff while even can be nonzero.
   For m', the odd yy response is nonzero. For 2', **an allowed even
   in-plane tensor nevertheless vanishes in the chosen model**;
   symmetry permission does NOT guarantee a nonzero physical response.
5. Wrong k-space TR action and replacing alpha/chi(-z) by a fixed-positive-
   eta negative-frequency retarded branch both fail deliberately.
6. Unit-tested analytic velocity from independent central momentum
   finite differences. Hexagonal BZ integration on an untransformed square
   k grid raises a clear error; C3/C6 pointwise tests are local smooth
   models and DO NOT claim a correct global periodic hexagonal crystal.

## Scope and unresolved source-vs-model questions

The 31/31 count is **generic-point local four-band Kubo symmetry covariance**.
It is NOT 31 complete symmetry-appropriate microscopic materials or 31
whole-Brillouin-zone experimental responses. The assumed source Table I
sector labeling and six prior conditional differences remain unresolved.
The current nondegenerate directed spectral response does not address
symmetry-enforced Kramers/PT doublets, finite-scattering conserving KME,
nonlocal magnetization current, full Wannier spin/position, or SI output.

21 geometry + 18 causal + 17 conventional spin + 24 vorticity/MPG
+ 20 three-group microscopic + 18 new tests = **118 tests**.
Prior accepted suite of 100 remains unchanged. No private AutoWann
files, licensed DFT input, publisher PDF, or material operator matrices
are included.

## Reproduce

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_group_microscopic.py

All counts and results are explicitly bounded model-method accomplishments.
See REPRODUCTION_PROGRESS.md for the three-track live status and blockers.
