# Round 24 | Kramers-degenerate equilibrium Zeeman–current Ward identity

**Scientific status:** independent bounded equilibrium-static method verification, NOT
the full original Cui et al. paper (still PARTIAL), NOT a finite-frequency or
dissipative magnetoelectric transport coefficient.

## Physics

In a fixed orthonormal Bloch basis with physical spin \(S_b\) independent of k,
and at finite temperature \(T>0\), consider
\(H_B(k)=H(k)-\kappa B_b S_b\) and \(J_a=\partial_{k_a}H(k)\).
The exact Fermi-matrix Fréchet derivative is

\[
D f_H[A] = \sum_{n,m} f[\epsilon_n,\epsilon_m]
|n\rangle\langle n|A|m\rangle\langle m|,
\]
\[
f[x,y]=\begin{cases}
\displaystyle\frac{f(x)-f(y)}{x-y},&x\ne y,\\
f'(x),&x=y,
\end{cases}
\qquad f(x)=\frac{1}{1+\exp[(x-\mu)/T]}.
\]

The complete *equilibrium static* field derivative of the current is

\[
\boxed{\chi_{JB}^{ab}(k;0)
=-\kappa\,\mathrm{Tr}[J_a\,D f_H[S_b]]
=-\kappa\,\partial_{k_a}\mathrm{Tr}[f(H(k))S_b].}
\]

This does **not** require a nondegenerate eigenbasis. At exact Kramers
degeneracy, the repeated-energy divided difference is \(f'\), including
**off-diagonal matrix elements within the degenerate subspace**.
Retaining only \(n=m\) Fermi-derivative entries would give a
gauge-dependent, physically incomplete static response.

For an exactly degenerate projector \(P_C\), the intra-cluster term is
\(-\kappa f'(\epsilon_C)\mathrm{Tr}[P_C J_a P_C S_b]\),
and the complementary inter-cluster term is the divided-difference sum
over \(C\ne D\). The split is gauge invariant under arbitrary \(U(d_C)\)
basis changes for exact degeneracies. Close but nondegenerate levels are
not flattened: the full divided-difference response remains continuous
without introducing a cluster approximation.

On a complete periodic reciprocal torus,
\(\int_{\mathrm{BZ}}\chi_{JB}^{ab}(k;0)\,d^2k=0\) for the assumed
k-independent \(S_b\); equilibrium cannot support a uniform bulk
charge current from a static spatially uniform Zeeman field.

## Numerical validations

- A genuinely TR-symmetric, triangular-periodic \(4\times4\)
  orbital × **physical-spin** Hamiltonian with two exact Kramers
  doublets at the zone center, and independent analytic
  Cartesian velocity \(\partial_kH\).
- Compare eigenspace Fréchet sum against a separate **spectral projector**
  contraction (avoids any degenerate single-band gauge).
- Independent uniform-B finite difference of equilibrium current,
  and momentum finite difference of equilibrium spin density.
- Arbitrary global-unitary covariance; explicit randomized within-doublet
  gauge, with a **negative control** showing the simple diagonal-only
  occupation derivative changes under a within-doublet rotation.
- Near-degeneracy continuity of the Fermi divided difference, symmetry
  under exchanging the energy pair, energy origin-shift invariance,
  reciprocal periodicity, TR at generic k, and fail-closed operator checks.
- Recovery of the earlier Round23 **nondegenerate** equilibrium Ward formula
  away from symmetry-enforced degeneracy. The old Round23 implementation
  continues to reject Kramers-degenerate states by design.

All numerical values are synthetic natural units. Independent tests are
required on both Python 3.11 and Python 3.12 before acceptance.

## Deliberate boundaries

No finite-frequency conserving-collision KME, no electromagnetic vector
potential/diamagnetic contacts, no realistic Wannier basis embedding or
k-dependent spin operators, no T=0 Fermi-surface distribution, no
non-Abelian dynamical Fermi-surface transport, and no SI response.
The extension controls an *equilibrium static limit* and provides an
essential negative oracle for subsequent non-equilibrium calculations.

## Reproduce

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_degenerate_static_ward.py

The original research paper's completion status remains **PARTIAL**.
