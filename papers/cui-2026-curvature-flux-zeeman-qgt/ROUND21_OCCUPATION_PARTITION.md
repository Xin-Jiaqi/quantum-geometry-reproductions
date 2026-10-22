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
