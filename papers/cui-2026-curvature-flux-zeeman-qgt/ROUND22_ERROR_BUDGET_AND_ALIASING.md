# Round 22 — BZ cancellation-conditioned error budget and finite-grid ensemble

**Scientific status: limited model-method; full Cui paper PARTIAL.**
**BZ numerical convergence remains UNVERIFIED** even if all tests pass.

## Why relative convergence of small net tensors is demanding

For the normalized BZ measure dµ and a FIXED complex observable tensor
component F(k), define I=∫F dµ, L=∫|F|dµ and retention r=|I|/L.

If an approximate numerical integrand satisfies the *unverified premise*
∫|δF|dµ ≤ eps*L, then by the triangle inequality

    |δI|/|I| ≤ eps/r,    whenever I!=0.

Thus to target 1% relative accuracy of I under this perturbation model,
the **necessary sufficient-for-this-bound** weighted-L1 tolerance is
eps ≤ 0.01*r. This is a *conditional mathematical budget*, NOT an estimate
of the actual quadrature error, a certified bound, or a real-material claim.

The Round20 fixed channels on the toy model's hexagonal order-5 mesh yield:

| Synthetic observable | r | Required eps for <=1% net bound |
| --- | ---: | ---: |
| 6, even yy | 0.0142484701 | 0.0001424847 = 0.014248% |
| 6m'm', even xy | 0.7887923373 | 0.0078879234 = 0.788792% |

Retentions here are themselves finite-quadrature approximations; the
conditional budgets illustrate **conditioning**, not achieved accuracy.

## Independent finite-grid comparison, NOT a confidence interval

The reusable kernel \`bz_error_budget.py\` implements:

1. Four independent shifts on a **coarse and refined** primitive reciprocal
   periodic midpoint mesh, using the actual four-band cluster Kubo.
2. Primitive tensor Gauss-Legendre integration and **hexagonal** first-BZ
   D6 symmetric Gauss integration using different quadrature domains.
3. Full complex chosen tensor component, never just a magnitude, and no
   magnetic-group post-projection of primitive estimates.
4. Fine-shift spread, coarse-to-fine midpoint mean difference, and
   finest-method spread, plus their maximum as an **empirical disagreement**.
5. Separate symmetry-allowed/forbidden classification and raw tensor
   projector residuals, without using symmetry residual as a convergence
   criterion.

**Crucial counterexample:** the exact primitive-cell average of
cos(8*q1) is zero, yet an UNshifted periodic 8×8 midpoint grid samples
the same Fourier phase at all nodes, returning a magnitude of one.
A quarter-index cell shift breaks this aliasing. Therefore finite-grid
agreement may be misleading and the method MUST NOT label numerical
convergence certified merely because its ensemble spread is small.

An empirical discrepancy can both under- and over-estimate the unknown
true integral error. It has **no asserted probabilistic coverage**, and
is **not** a rigorous integration-error upper bound. The output always
reports "UNVERIFIED" for physical convergence.

## Scope

- Model: same periodic synthetic spin × orbital 4-band Hamiltonian as
  Round17–20; same chemical potential, temperature and complex drive.
- Original Cui et al. complete Table I, experimental KME, SI units,
  conserving metallic collisions, real Wannier spin/position and full
  paper reproduction remain outside this numerical method.
- Do not introduce a new material or VASP calculation at this gate.

## Reproduce

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_round22.py
    python scripts/audit_round22.py --deep

11 new tests + 167 current baseline = **178 target**, pending remote CI. The
frozen JSON will record exact physical channels, mesh resolutions,
the observed empirical complex disagreement, conditional L1 budgets,
and the explicit non-certification.



## Preliminary independent model calculation (pending integrated branch CI)

This method was first tested in isolation against the Round20 baseline
(170/170 tests successful on both Python versions in Actions run
39748473393). While that evaluation was underway, separately verified
**Round21 occupation_partition.py** entered main (eight additional tests,
167/167 PASS). This error-budget capability is therefore integrated as
**Round22**, leaving Round21 and its frozen outputs untouched.

For group 6 even yy, the hexagon order-5 magnitude is 4.96237e-5,
whereas the finest 9x9 shifted midpoint estimate is around 5e-7.
Fine-grid shift spread is 2.37929e-7; the largest independent-method
complex disagreement is 4.90624e-5, approximately 98.87% of the finite
hex5 anchor magnitude.

For group 6m'm' even xy, the hex5 magnitude is 8.39843e-7.
Fine-grid shifted midpoint spread is 4.43700e-7 and the empirical
disagreement relative to finite hex5 is 52.83%.

These are FINITE-grid differences, NOT rigorous errors relative to the
unknown continuum integral. Exact full-tensor agreement with a separate
independent spinor spectral Kubo implementation was observed for sampled
grids, not for the continuum limit.

The integrated Round22 suite target is **178 tests** on Python 3.11
and 3.12; record final CI results only after observing them.
