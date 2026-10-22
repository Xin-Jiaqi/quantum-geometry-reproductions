# Round 23 | Equilibrium static Zeeman-current Ward identity

**Scope:** exact clean, static, fixed orthonormal Bloch basis with a k-independent
physical-spin operator. The original Cui et al. paper stays **PARTIAL**.
This is an independently derived *physical consistency control*, not a
reproduction of an additional author equation or complete finite-frequency
KME.

## Why the prior band-only response can be misleading

Let H_B(k)=H(k)-kappa B_b S_b and J_a(k)=partial H(k)/partial k_a.
Assume S_b is independent of k in the declared basis. At finite T>0
and fixed chemical potential the equilibrium Bloch current is
j_a(B)=integral_BZ Tr[f(H_B(k)) J_a(k)].
As J_a=partial_k H_B and f is the derivative of a scalar thermodynamic
primitive F, we have Tr[f(H_B)J_a]=partial_{k_a}Tr[F(H_B)].
Its integral on the reciprocal torus is zero for any uniform B.

Differentiating at B=0 gives the *local* static Ward identity

 chi_JB,static^(ab)(k)
   = chi_interband^(ab)(k)+chi_occupation^(ab)(k)
   = -kappa partial_{k_a} Tr[f(H(k)) S_b],

with
 chi_interband = -kappa sum_(n!=m)
 [(f_n-f_m)/(eps_n-eps_m)] v^a_nm S^b_mn
and
 chi_occupation = -kappa sum_n f'(eps_n) v^a_nn S^b_nn.

The second term is the occupation derivative / intraband static piece;
it is *not* a phenomenological damping term. It must not be omitted
when promoting a bare interband Kubo response to a complete
equilibrium static current response in a metal at finite T.
It vanishes exponentially if mu lies well inside a gap and T→0.

## Independent checks

1. Direct spectral response from the above two-term formula and the
   existing directed interband Kubo at z=0 agree for nondegenerate bands.
2. Finite-Zeeman-field central difference of Tr[f(H-kappa B S)J]
   agrees with the full interband+diagonal static response.
3. Independent finite-momentum difference of the equilibrium spin density
   reproduces minus kappa times the complete response.
4. Constant orbital-unitary covariance, realistic nonzero f', strongly
   gapped low-T limit, and fail-closed degeneracy tests.
5. Full periodic BZ is a torus. A periodic central finite difference of
   spin density telescopes to exactly zero on *every* complete midpoint
   grid; a sampled spectral derivative need not do so at finite order.
   Any nonzero finite-grid static current is numerical integration error,
   not an equilibrium KME observable.

Run:

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_equilibrium_ward.py

**Exclusions:** arbitrary k-dependent physical spin/Wannier embedding,
contact corrections, finite-q electromagnetic gauge, inelastic transport,
finite-frequency driven metallic KME, SI response, and genuine
degenerate-band non-Abelian response in this *static spectral split*.
Degenerate models are explicitly rejected rather than assigned
gauge-dependent diagonal Fermi derivatives.

The equilibrium no-current theorem places a strong negative control
on any static B-induced current in a periodic bulk crystal.
It does not prohibit finite-frequency response, non-equilibrium transport,
edge currents, magnetization circulating currents, or spin currents.


## Actual dual-version CI and frozen static equilibrium data

GitHub Actions run **39756688525**, Python 3.11 and Python 3.12, both
executed **188/188 tests, PASS** and the independent static audit.
Python 3.12's raw JSON is frozen at
results/frozen/round23_equilibrium_ward_audit.json.

At k=(0.47,-0.31) in the full four-band physical-spin/orbital synthetic
model with mu=2.3, T=0.3 (natural units):

- Max residual, full analytic static susceptibility vs independent
  finite uniform Zeeman field derivative: **7.880842838e-11**.
- Max residual, analytic static susceptibility vs independent finite
  momentum derivative of equilibrium spin density: **1.253953885e-10**.
- Interband and occupation-derivative terms are both substantial; for
  component (a=y,b=x), interband = +0.0140794943 and occupation
  derivative = +0.2260539938, giving full +0.2401334881.
  Thus static band-only Kubo would miss the dominant diagonal term.

| Primitive midpoint k mesh | Max full static spectral BZ mean | Exact discrete periodic telescope |
|---|---:|---:|
| 4 x 4 | 1.157326428e-2 | 4.337e-19 |
| 8 x 8 | 7.062745359e-5 | 2.602e-18 |
| 12 x 12 | 8.685197330e-7 | 7.999e-18 |

The spectral integral approaches the exact periodic no-current identity
as the mesh is refined. This finite sequence is a convergence *illustration*,
not a rigorous error bound. The exactly vanishing discrete telescope is a
separate (and mathematically exact) finite-grid summation-by-parts oracle:
it must NOT be substituted as a value of the spectral Kubo quadrature.

**Physical conclusion:** a uniform static B perturbation cannot produce
an equilibrium bulk charge current on the complete periodic BZ under the
stated assumptions. An apparent nonzero small static interband value is
not automatically an equilibrium magnetoelectric transport effect.
The dynamic and dissipative KME remain open.
