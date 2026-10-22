# Round 20 | Physical origins of slow BZ response convergence

**Scientific state:** model-method under review; full paper PARTIAL. All
numbers below are natural-unit *synthetic four-band orbital × physical-spin*
model values at z=0.071+0.004i, mu=0.13, T=0.17, from an independently
reconstructed implementation, not real-material or author data.

## Pair-resolved Kubo partition

Let chi_even(k) = [chi_JB(k,z)+chi_JB(k,-z)]/2. The directed inter-cluster
transition sum has an exact two-way split:

- cross-mu: the pair energies lie on opposite sides of mu;
- same-side: the two energies lie on the same side of mu.

These two terms reconstruct the original Kubo tensor exactly at each k.
The near-resonant (detuning <= 2 eta) and small-gap (gap <= 0.02)
subsets are OVERLAPPING diagnostics, not extra independent channels.

The cancellation ratio for a fixed Cartesian tensor component ab is
R_ab = sum_i w_i |chi_even^ab(k_i)| / |sum_i w_i chi_even^ab(k_i)|.
Large R implies severe sign/complex-phase cancellation, not a certified
numerical error bound. A symmetry zero can make R diverge; report only
a physically allowed component when drawing conclusions.

## Independently computed results: true first hexagonal BZ, order 12

**Sixfold group 6, allowed xx component:**

- net chi_even(xx) = 4.83881e-7 + 3.99709e-8 i;
- integral of absolute pointwise kernel = 3.36024e-3;
- cancellation ratio about 6921.
- Cross-mu absolute integral = 3.36023e-3 but its integrated net is
  only 3.74927e-9 + 3.61778e-10 i.
- Same-side absolute integral = 1.12897e-6, but its net is
  4.80132e-7 + 3.96091e-8 i.
- The largest Fermi occupation difference among sampled pairs with
  energy pole detuning <= 2*eta is ~1.08e-6. A small energy denominator
  does NOT imply a large response if the initial/final occupations match.

For this finite-grid toy model the dominant local inter-manifold virtual
terms cancel by thousands-fold, leaving a small thermally weighted
same-side contribution. This is a clue, NOT proof that the converged
physical response has that decomposition.

**Group 6m'm', allowed xy component:**

- net chi_even(xy) = 2.51817e-7 + 6.77616e-7 i;
- absolute pointwise kernel average = 1.04219e-6;
- cancellation ratio only about 1.44;
- cross-mu net is ~8.0e-10; same-side net is
  2.51016e-7 + 6.77615e-7 i.

Thus different symmetry sectors have very different cancellation behavior.
Do not claim all BZ convergence problems are simply small gaps or giant
phase cancellation.

## Verification contract and exclusions

1. Reconstruct positive/negative full complex frequency Kubo tensors
   at fixed k from individual pair contributions, independently comparing
   against the pre-existing cluster_interband_kubo.
2. Verify sum(cross_mu,same_side)==total and do NOT double count
   resonance/small-gap diagnostic subsets.
3. Integrate without post-projecting into allowed magnetic tensor space.
   Check against the previous hexagonal and primitive BZ modules at the
   SAME k-grid and weights.
4. Fix H(k), mu, z while scanning temperature; report band gaps unchanged
   and the relative occupation contributions.
5. Explicitly distinguish net signed complex response, average magnitude,
   min gap, frequency detuning and occupation difference.

Open: true high-accuracy BZ convergence, an analytic occupied/full-band
sum rule, non-Abelian intraband metallic KME, collisions and gauge contacts,
validated Wannier physical-spin operators and source Table I/author Eq.(19)
normalization. No SI or real material result is asserted.
