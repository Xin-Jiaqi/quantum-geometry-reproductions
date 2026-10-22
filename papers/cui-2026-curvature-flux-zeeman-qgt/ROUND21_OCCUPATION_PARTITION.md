# Round 21 | Occupation-sector resolution of the Cui Flux k-space response

**Status:** bounded synthetic-model diagnostic; whole-paper PARTIAL.
This is deliberately a *complement* to the independently verified
Round20 k-resolved hotspot module (bz_hotspots.py, 159 public tests).
It does NOT replace that module or rerun all its hotspot and symmetry
audits. The new question is which band-pair OCCUPATION sector survives
after k-integration.

## Exact pair partition

For each k, every directed inter-cluster clean Kubo term is assigned to
exactly one of two disjoint categories:

- **Cross-mu:** (eps_n - mu)(eps_m - mu) < 0.
- **Same-side:** both band energies on the same side of mu.

Thus chi_even(k) = chi_cross_mu(k)+chi_same_side(k) as a strict algebraic
identity. The same-side population difference may be exponentially small,
but one must not discard it merely because the local contribution appears
negligible relative to other virtual interband transitions.

Near-resonant and small-gap terms are OVERLAPPING annotations of the same
pairs; they are *not* added again to the complete cross/same partition.
The near-pole criterion uses the full complex drive z, and occupation
differences are reported separately from denominator magnitude.

## Independently reconstructed toy-model evidence

At z=0.071+0.004i, mu=0.13, T=0.17, synthetic 4-band real-spin /
orbital model, hexagonal first BZ D6-symmetrized order12:

**Group 6, even xx (allowed):**
- Net 4.83881e-7+3.99709e-8 i.
- Weighted absolute local-k kernel = 3.36024e-3; cancellation ratio 6921.
- Cross-mu weighted absolute kernel = 3.36023e-3, but signed net
  ~3.74927e-9+3.61778e-10 i.
- Same-side weighted absolute kernel ~1.12897e-6, but signed net
  ~4.80132e-7+3.96091e-8 i.

**Group 6m'm', even xy (allowed):**
- Net 2.51817e-7+6.77616e-7 i.
- Weighted absolute kernel 1.04219e-6, ratio 1.44.
- Cross-mu net ~8.00e-10, while the same-side net is
  ~2.51016e-7+6.77615e-7 i.

Therefore the microscopic origins of weak net response differ across
allowed magnetic sectors. The tiny local occupation difference of pairs
near Re z~0.071 can be ~1e-6; small denominators alone do not establish
an important contribution.

## Validation and interpretation

- The new pair sum separately reconstructs existing complete chi(z)
  and chi(-z) from cluster_interband_kubo, then exact even decomposition.
- The cross-mu plus same-side sectors reconstruct the full tensor
  pointwise and after independent hexagonal/primitive BZ sampling.
- Fixed-H temperature scans isolate the Fermi-factor change while
  energy gaps and velocity operators remain unchanged.
- The existing bz_hotspots.py is the authoritative spatial hotspot
  sampler; this add-on ONLY publishes thermal occupation sectors.
- Test count: 8 additional tests on top of the accepted 159-test
  public suite, to be certified by both Python 3.11 and 3.12 CI.

**Strict caveat:** these are finite-order toy-model diagnostics. The
cross-mu near cancellation is not yet proven as an exact sum rule,
and the remaining same-side contribution is NOT a converged material
current. No Fermi-surface Drude term, conserving collision, real Wannier
physical spin/position, SI response, or whole-paper PASS is asserted.

## Reproduce

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_occupation_partition.py

This package must be accepted as **bounded model method**, not as
replacement of the already accepted Round20 module.


## Executed GitHub CI and frozen numerical reference (2026-10-22)

Actions run **39748027692**, Python 3.11 and Python 3.12:
**167/167 tests PASS**, and both completed the independent scripted
sector audit. The actual Python 3.12 audit JSON (not a manually
constructed model benchmark) is frozen at
results/frozen/round21_fermi_sector_audit.json.

- Group 6, hexagon order8, allowed even xx: net
  **4.017570968e-7 + 2.475013118e-9 i**, cancellation ratio 8265.23.
  Its cross-mu net is **-2.007035130e-7+3.539161094e-10 i**,
  while its same-side net is **6.024606098e-7+2.121097009e-9 i**.
  Their sum exactly reconstructs the full finite-grid response.
- Group 6m'm', hexagon order8, allowed even xy: net
  **3.980107636e-7+4.848696899e-7 i**; cancellation ratio 1.4903.
  Cross-mu net ~8.00051e-10; same-side net
  **3.972107126e-7+4.848696706e-7 i**.
- Group 6, primitive Gauss order12, even xx: net
  **2.256915657e-7+3.547745026e-7 i**; cancellation ratio 7959.97.
  Its full tensor differs from hex-order8: **no BZ convergence claim**.
- At fixed k=(0.5145,-1.0456), same-side xx changes from
  ~1.11296e-13 at T=.08 to 1.50417e-7 at T=.17 and 5.92000e-6
  at T=.25, whereas the minimum intercluster gap stays
  **0.17474288932364512**. This isolates the thermal occupation
  factor without changing H or velocity.

These findings are independent of and *complement* the prior
159-test Round20 hotspot acceptance. The only newly promoted method
is the exact cross-mu versus same-side thermal pair partition, NOT a
second competing hotspot implementation or a converged experiment.
