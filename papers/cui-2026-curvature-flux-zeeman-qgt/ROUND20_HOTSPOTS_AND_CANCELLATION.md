# Round 20: k-resolved spectral gaps, occupied transitions and BZ cancellation

**Scientific acceptance:** synthetic nondegenerate/cluster interband diagnostics VERIFIED;
**numerical BZ convergence remains UNVERIFIED; full-paper status PARTIAL**.

## Independent physical formula

For a directed inter-cluster transition \(\Delta=\epsilon_m-\epsilon_n\)
and \(A_{nm}^{ab}=(f_n-f_m)v^a_{nm}s^b_{mn}\), the clean complex-frequency
electric-current / Zeeman-field response obeys

\[
\chi^\mathrm{even}_{ab}(k,z)=-\kappa\sum_{n,m\in\mathrm{different~clusters}}
\frac{A_{nm}^{ab}\Delta}{z^2-\Delta^2},\qquad
\chi^\mathrm{odd}_{ab}(k,z)=-\kappa\sum_{n,m\in\mathrm{different~clusters}}
\frac{A_{nm}^{ab}z}{z^2-\Delta^2}.
\]

This is an algebraically independent numerical evaluator compared with
separately computing \(\chi(k,z)\) and \(\chi(k,-z)\) through the prior
Kubo implementation. Both agree at arbitrary k and after integration.

**Important corrected error:** an early implementation wrongly used a
reporting threshold on \(|f_n-f_m|\) to delete small-occupation physical
terms. That failed independent Kubo regression. The fixed implementation
includes every nonzero occupation difference in the physical tensor.
Thresholds affect only **reported active-transition counts/detunings**;
a new negative-control test ensures changing the reporting cutoff cannot
change the tensor.

## Cancellation and hotspot definition

For a *fixed Cartesian tensor component*, define
\(N_{ab}=\sum_k w_k\chi_{ab}(k)\),
\(L_{ab}=\sum_k w_k|\chi_{ab}(k)|\),
and cancellation \(C_{ab}=1-|N_{ab}|/L_{ab}\).
Report top-5%-node share of L by ranking **weighted** local absolute
contributions, not raw integrand values alone. Full complex tensors are
always preserved without any symmetry projection.

Symmetry-equivalent components can have equal integrated norms. A first
audit automatically picked xy or yx depending on float roundoff and
misleadingly produced different hotspot statistics under Python 3.11 vs
3.12. This was caught and fixed by explicitly specifying the physical
observable (group 6: yy; group 6m'm': xy) and using a deterministic
relative-tie rule for optional auto-selection.

## Frozen illustrative model results

Drive \(z=0.071+0.004i\), chemical potential \(\mu=0.13\),
temperature \(T=0.17\), natural units.

| Magnetic group / numerical grid | Observable | Integrated local magnitude L | Retention 1-C | Cancellation C | Top 5% weighted share |
| --- | --- | ---: | ---: | ---: | ---: |
| 6, hexagonal Gauss order 5 | yy even | 3.48274e-3 | 0.0142485 | **98.575%** | 23.12% |
| 6, primitive Gauss order 10 | yy even | 3.07213e-3 | 0.00031818 | 99.968% | 25.85% |
| 6, primitive shifted midpoint 12×12 | yy even | 3.26254e-3 | 0.00018323 | 99.982% | 18.57% |
| 6m'm', hexagonal Gauss order 5 | xy even | 1.06472e-6 | 0.788792 | **21.12%** | 31.22% |

The net group-6 yy component (~4.96237e-5 at hexagonal order 5) is
much smaller than its absolute local response integral (~3.48e-3).
Small absolute quadrature errors therefore become large **relative**
errors in the cancellation-dominated net. The two groups show
qualitatively different cancellation fractions; the mechanism cannot
be generalized from one magnetic point group to another.

Minimum sampled intercluster gaps of order 0.006–0.011 coexist with
minimum *active-transition detunings* around 4.4 in these grids.
Here 'active' only describes pairs passing the diagnostic threshold
\(|f_n-f_m|>10^{-5}\), **not** a cutoff in the actual Kubo tensor.
The small minimum gap alone does not demonstrate an occupied optical
resonance at Re z=0.071; weakly occupied near-degenerate transitions
are still retained and require separate study.

A hexagonal-order-4 group-6 T scan changes the maximum even response
from 3.43540e-4 (T=0.10) to 3.42995e-4 (T=0.17) and 2.54306e-4
(T=0.35). These are **finite-grid occupation sensitivity checks**, not
converged thermodynamic response coefficients.

## Numerical validation and reproducibility

- **14 new tests; 159 total**. Python 3.11/3.12 both PASS at PR head
  f15b5e542c19c0c477499c35ba136d4fd3885329.
- Max cross-Python complex tensor element difference: 3.62e-18;
  max selected scalar-diagnostic difference: 7.94e-13.
- Pointwise independently regrouped pair denominator vs prior Kubo,
  hexagonal integral, unprojected primitive integral, global energy
  shift with chemical potential, threshold isolation and stable symmetry
  channel selection are tested.
- Full JSON audit, node hotspots, gaps, occupations, full even/odd complex
  tensors, and cross-Python provenance are stored at
  results/frozen/round20_khotspot_cancellation.json.

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_round20.py
    python scripts/audit_round20.py --deep

## Unresolved questions

High-order BZ convergence of the small net tensor is **not established**.
A sampling minimum gap is not a global band minimum. The diagnostic's
active-transition detuning does not prove all weakly occupied
near-degenerate contributions irrelevant. No metallic intraband response,
conserving collision, true Wannier spin/position matrices, real material,
SI response, or full source Table I/Eq.(19) is verified by this work.
