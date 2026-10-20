# BSF × NMEE: stage 3 scientific acceptance

## 1. New result within a fixed synthetic model

We kept the *microscopic real-space displacement-generated hopping recipe*
from Stage 2 and varied only relative sliding. Parameterize the fractional
shift by `tau(u) = u(0.23,0.41)` and fix `t_perp = 0.020 eV`, opposite layer
Rashba SOC `0.020 eV`, intrinsic SOC `0.010 eV`, `mu=0.02541381265 eV`,
zero layer potential difference, T=0, and a Zeeman finite-difference step of
0.25 T. The response is `Gamma_{xz;x}` in `(nm V T)^-1` under the
**frozen-occupation, interband static susceptibility** protocol, **not** the
unresolved driven-DC Eq. (1) of Jia et al., PRL 136, 226302 (2026).

A 16-point scan reveals a sign reversal without a sampled insulating gap
closure; selected results are:

| Path u | Static Gamma_xz;x | Sampled minimum distance of mu to the bands |
|---:|---:|---:|
| 0 | approximately 0 | 2.10 meV |
| 0.10 | -1.6981e-5 | 2.43 meV |
| 0.20 | +9.32e-7 | 3.37 meV |
| 0.40 | +9.8502e-5 | 6.40 meV |
| 1.00 | +1.8289e-4 | 14.47 meV |

The zero crossing is bracketed by **u=0.1973877–0.1974121** at the fixed
coarse quadrature (`u≈0.197400`); it is a **non-symmetry-protected, model-
dependent cancellation** rather than a topological phase transition.

For `t_perp = 0.028 eV`, the original fixed chemical potential falls outside
the insulating gap at `u=0`. **This is why the path scan uses 0.020 eV.** The
Stage 2 single-shift 0.028 eV benchmark is not being revised or compared as a
same-parameter point.

## 2. Integrator and independent convergence

The BZ is partitioned into two non-overlapping K/K′ disks (radius 1.5 /nm)
plus the exterior. The two resolutions (exterior 41² + valley 35×64,
exterior 61² + valley 60×80) produce, for `u=0.10,0.20,0.40,1.00`,
absolute amplitude differences below **3e-9**. Each result also remains of
the same sign with the `0.125 T` field step (tested u=0.10,0.20,0.40).
The low-resolution exterior measure differs from the exact BZ area by
approximately 0.15%; the exterior response itself is negligible here and
the fully recomputed denser quadrature confirms the reported signal.

This is a **numerical quadrature acceptance** on sampled models, not a proof
of global bandgap throughout the continuous toroidal BZ; the absolute
thermodynamic/experimental magnitude remains unverified.

## 3. Valley-resolved and orbital-resolved mechanisms

At u=0.10, each valley contributes approximately `-8.49e-6` and at u=0.40
each contributes approximately `+4.925e-5`. The rest of the BZ is of order 10^-9 or less.
**Consequently the sign reversal is not K/K′ valley cancellation:** both
valleys reverse together. Radial shell checks localize the main reversal
within approximately 0.25 /nm of each valley center. These patch
contributions depend on the stated valley partition and are not independent
material observables.

We separately turn off A↔B *cross-sublattice* interlayer hopping. Because
removing channels also changes total tunneling strength, a stronger control
matches the BZ-averaged `Tr[T†(k)T(k)]` (41² k mesh) by rescaling the remaining
same-sublattice hopping:

| u | full hopping | same-only, equal BZ-averaged hopping norm |
|---:|---:|---:|
| 0.40 | +9.850e-5 | +2.57e-8 |
| 1.00 | +1.829e-4 | +3.34e-8 |

Thus **orbital mixing matters strongly even after removing the trivial
change in the hopping norm**. However, norm matching does not fix the full
band dispersion, orbital-resolved spectral weights or local gaps: the result
is a strong *model mechanism clue*, not a universal selection rule, nor proof
that a uniquely defined Zeeman-metric contribution is responsible.

## 4. Symmetry and physics boundaries

The same finite-range hopping function generates all shifts independently.
The physical horizontal-mirror/layer-exchange relation gives
`Gamma_total(-tau)=-Gamma_total(tau)` and
`Gamma_top(-tau)=-Gamma_bottom(tau)`, numerically to around 1e-16.
This symmetry **does not force** `Gamma(tau(u))` to keep one sign between the
mirror-related endpoints. Conversely, it enforces `Gamma(0)=0` provided the
high-symmetry structure is insulating.

**NOT established:** physical BSF minima or energy barriers; existence of a
ferroelectric branch; generality across chemically distinct models;
impurity/non-equilibrium DC responses; a gauge-invariant band/orbital
attribution; actual Wannier materials, layer position operator or measured
amplitudes. A sign change by itself is **not PRL-level novelty**. Next
scientific step is to identify a parameter-robust geometrical or symmetry
principle and eliminate mere spectral redistribution explanations.

Source: `results/sliding_stage3_sign_reversal.json` (all parameters, raw
layer-resolved values, full k integration tests, gap diagnostics),
`scripts/run_stage3_sign_reversal.py` (recompute), and
`tests/test_sliding_sign_reversal.py` (test gates).
