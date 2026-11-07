# Bilayer static EB spin response: polarization, powers and coupling

Date: 2026-11-03. [Short Chinese explanation](BILAYER_KM_FINDINGS_ZH.md).
This audit supplies falsifiable model results, not a material prediction or a
claim that the complete driven-DC NMEE problem has been solved.

## Observable and model contract

Use the insulating, zero-temperature, static branch throughout. For i,j=x,y,

\[
\delta s_i=C_{ij}E_jB_z,\qquad C_{ij}=\partial_{B_z}\chi_{ij}|_{B_z=0},
\quad
\chi_{ij}=-2\operatorname{Im}\int[dk]\sum_{v,c}
\frac{(\sigma_i)_{vc}(\partial_{k_j}H)_{cv}}{(\epsilon_v-\epsilon_c)^2}.
\]

Here [dk]=d²k/(2π)², v spans four occupied bands and c the four empty bands.
The observable is total spin of the bilayer, not layer-resolved spin or spin
current. This uses the repository convention for electron charge and Pauli spin;
C is reported in legacy 1/(nm V T) units. A physical spin angular momentum density
requires the factor ħ/2. The Zeeman convention is H(B)=H(0)−g μB B σz with g=10;
this effective g includes any factor convention and must not be silently replaced
by the usual −g_phys μB B σz/2. No orbital magnetic coupling is included.

This definition agrees numerically with `symmetry_scaling.static_inplane_matrix`.
It does not identify the static coefficient with the printed PRL intrinsic
driven-DC coefficient. Scattering, finite-temperature occupations and order of
frequency/relaxation limits are not addressed. Do not mix those branches.

Default monolayer: t=0.85 eV, bond=0.142 nm, sublattice potential=10 meV,
intrinsic SOC=10 meV, Rashba SOC=20 meV. The opposite layer is σz Hmono σz,
which reverses Rashba. Interlayer hopping is the existing displacement-dependent
five-shell Gaussian model (see stored `parameters` for all values), t_perp=5 meV.
The chemical potential is fixed inside the gap. Fractional displacement is
d_frac=u(0.23,0.41); u is amplitude, not angle. Angle controls rotate d in Cartesian
space before converting back to fractional coordinates.

The nonzero staggered and Rashba terms invalidate the full graphene/LG80 premise.
Off-grid Hamiltonian tests verify monolayer C3 and a vertical mirror, and show that
its horizontal mirror is broken. The bilayer reference u=0 restores horizontal
mirror by exchanging the opposite-Rashba layers. At generic u the mirror relates
two different structures ±u; C3 likewise relates rotated displacements, not a
symmetry of each arbitrary displaced structure. Time reversal holds at B=0.
Residuals: C3 <4e−16 eV, layer-mirror relation <2e−18 eV, TR zero at tested k.
The ordinary graphene inversion test fails by 20 meV for the default Hamiltonian.

The direct-stacking control uses identical full monolayer Hamiltonians, and has
a large nonzero response already at u=0; raw powers there cannot be interpreted
as onset powers without subtracting the reference. A separate direct-stacked
graphene control with staggered=Rashba=0 preserves inversion at any tested slide
and gives zero dipole and net EB response to numerical precision.

## Added polarization observable

For fixed point orbitals at z=±d/2 and neutral, equal rigid ionic backgrounds,

\[
\Delta N=A_{cell}\int[dk]\operatorname{Tr}[P_{occ}(P_{top}-P_{bottom})],
\quad D_z/(ed)=-\Delta N/2,\quad P_z^{2D}=D_z/A_{cell}.
\]

The minus sign is the electron charge, e>0. Since d is not specified in the
Hamiltonian, store D_z/(ed) per cell and P_z/(ed) per nm²; do not invent a layer
spacing or normalize by vacuum thickness. This is a well-defined electronic
out-of-plane dipole in the point-orbital model. It is not an in-plane Berry-phase
polarization, a relaxed-ion DFT polarization, or evidence of ferroelectricity.
A finite-difference derivative of occupied energy with respect to layer bias
independently reproduces ΔN (Hellmann–Feynman test).

## Result 1: local powers and their domain of validity

Define T=Cxx+Cyy and D=Cxx−Cyy. The symmetric traceless sector is the pair
(D,Cxy+Cyx), so one component D may vanish for a special sliding direction.
These spatial tensor combinations are unrelated to sums/differences between
opposite magnetic domains; this model has no magnetic order parameter.

For u=0.025→0.05, n=log2|Q(2u)/Q(u)| gives:

| Quantity | Local exponent | Interpretation |
|---|---:|---|
| Electronic Pz | 2.997 | cubic onset |
| T | 2.996 | cubic onset |
| D | 0.992 | linear onset |

The analytic selection rule is local about u=0 and assumes an open gap. The
horizontal mirror sends u→−u and Pz,C→−Pz,−C, so only odd powers occur. C3 makes
Pz and T in-plane scalars, excluding a linear invariant; the first permitted
odd polynomial is cubic. The rank-two traceless sector can couple linearly to
the displacement under C3 (angular momenta 2 and −1 coincide modulo 3).
This establishes allowed lowest orders, not nonzero coefficients or universality
at large slide. Symmetry constrains the same interlayer mixing that generates
the microscopic coefficients; they are not independent competing explanations.

Rashba, intrinsic SOC and staggered potential separately varied by ±20% retain
gaps and give scalar exponents 2.993–2.997 and D exponents 0.980–0.994. Six rotated
paths at 0,30,60,90,120,180 degrees retain the scalar cubic trend and directional
linear trend at the sampled amplitudes; this finite set does not exclude special
zero directions. Changing the first intralayer bond by 8 meV breaks C3 while
preserving the mirror relation. At u=0.003125→0.00625 the exponents become 0.995
(Pz), 0.993 (T), 1.000 (D). At larger u, competing linear and cubic terms cause
zeros and misleading noninteger two-point slopes. This is why the extra small-u
control was necessary.

At u=0.025 the electronic D_z/(ed) is only −1.12e−11 per cell; at u=0.1 it is
−7.09e−10. The normalized figure establishes powers, not experimental visibility
or an enhancement factor comparing quantities with different units.

## Result 2: what interlayer hopping does

At fixed u=0.2, turning off all interlayer hopping gives C≈0. From t_perp=1→2 meV,
D changes by a factor 4.006, consistent with a leading quadratic coupling effect.
In the fixed local A/B orbital basis, compare hopping between equal sublattices
(AA,BB) and unequal sublattices (AB,BA):

| Hopping model at 1 meV | D [1/(nm V T)] |
|---|---:|
| Both paths | −1.11815e−7 |
| Equal only | +3.37585e−11 |
| Unequal only | +2.46060e−10 |
| Both, relative sign reversed | +1.12508e−7 |

The full response is not the sum of the separately switched-on models. The large
directional difference is consistent with interference between the two paths
at weak coupling. Reversing both signs together is a layer gauge change and
leaves observables invariant (tested); reversing only one is a physically
different Hamiltonian. These controls establish a mechanism in this chosen
orbital model, not a unique gauge-independent decomposition of a single
observable into named quantum geometries. Small isolated-channel values are not
used as robust physical claims.

For t_perp=5 meV, Cxx changes from −7.97e−7 at u=0.2 to +2.17e−6 at u=0.4.
Cxx=(T+D)/2: a growing positive isotropic part overtakes the negative directional
part. Both Pz and the gap retain their sign/open status; this is neither required
polarization reversal nor a demonstrated topological transition.

## Result 3: energy changes alone cannot explain the response

Construct a counterfactual Hflat=Δ(I/2−Pocc), Δ=30 meV, for every k and u. Every
spectrum is exactly ±15 meV, while the zero-field occupied projector, layer
dipole and physical spin operators are retained. The velocity is computed
consistently as −Δ ∂k Pocc, not borrowed from the original Hamiltonian. Apply
the Zeeman perturbation after flattening. Tests check identical projectors,
fixed spectra and analytic velocities against finite differences.

The flattened C still varies with u, with cubic T and linear D at small u.
Therefore changes of electronic states suffice for sliding-dependent response
even when energies are identical. However, Cxx remains negative through u=0.4
for the flattened 5 and 20 meV coupling cases, while the physical cases change
sign. The observed physical sign reversal involves energy weights as well as
state geometry. Flattening is generally a long-range, artificial Hamiltonian;
it is neither an experimental control nor an additive “pure geometric” part
of the original C. The comparison does not isolate a unique named QGT.

## Numerical evidence and reproducibility

Results: `bilayer_stage8_audit.json` (main scan), `bilayer_stage8_controls.json`
(small-u broken-C3, direct reference, flattened convergence). Both in `results/`.
All cases are deterministic, without stochastic sampling or statistical error bars.

The valley/exterior quadrature is closed under C3 by independently evaluating
k,Rk,R²k with one-third weights. No observable is symmetrized or fitted to an
assumed power. This also integrates the explicitly C3-broken Hamiltonian.
The original exterior grid can leak spurious linear terms into tiny cubic
scalars; preserving the integration domain symmetry removes that leakage.
Grid refinement remains necessary: closure is not a replacement for convergence.

Base grid (61,60,96) has 45,513 evaluated points; fine (101,100,144) has 116,427.
At u=.025,.1,.4, relative base/fine changes are <0.062% for Pz, <0.00024% for T,
and <0.000017% for D. Halving B from .125 T to .0625 T changes T,D by <0.00081%.
For broken C3 at the smallest u, refinement changes Pz,T,D by <0.070%.
For the flattened controls at u=.1 and both couplings, grid and B-step changes
are <0.0053%. The raw quadrature area errors are retained in JSON metadata;
no hidden renormalization is applied. The default sliding scan's sampled
indirect gap stays above 28.78 meV. These are numerical checks, not rigorous
global analytic gap bounds.

From repository root, after installing the dependencies in the study README:

```sh
export PYTHONPATH="$PWD/papers/jia-2026-nonlin-magnetoelectric-edelstein:$PWD/studies/bsf-nmee-connection"
export OPENBLAS_NUM_THREADS=1
python3 studies/bsf-nmee-connection/scripts/run_stage8_audit.py
python3 studies/bsf-nmee-connection/scripts/run_stage8_controls.py
python3 studies/bsf-nmee-connection/scripts/plot_stage8_audit.py
python3 -m unittest discover -s studies/bsf-nmee-connection/tests -v
```

The explicit PYTHONPATH prevents a different editable installation from silently
supplying the physics kernel. `--quick` writes a separate quick JSON, which must
not be used as production evidence. Old Stage 1–7 manifests remain untouched.
The Stage 8 manifest freezes sources, data and figure exports; refresh it only
after the corresponding calculations/checks, never to suppress a failing test.

Local validation: 63/63 study tests pass, with the actual pinned v0.18.0 symmetry
toolkit and no skipped tests. The initial local run lacked that toolkit and
exposed the historical README hash protection; the historical README was restored,
and the new entry point was placed in `STAGE8_README.md` and the repository root.
No historical manifest was rewritten. Figure source preflight has zero failures,
and PNG rendering was inspected. Code review checked derivative/charge conventions,
degenerate-subspace projector treatment, legacy-kernel agreement, manifest scope,
bounded batch memory and the absence of private machine paths or credentials.

## Novelty boundary and next decision

This supports a candidate question: which structural symmetry conditions allow
an EB spin readout linear in a displacement whose charge dipole is only cubic?
The power distinction and its explicit failure on breaking C3 are useful
evidence. Alone they remain a selection rule plus one model, not a demonstrated
PRL-level new mechanism or an 80-layer-group classification.

Relevant prior art, with scope kept distinct:

- [Jia et al., NMEE](https://arxiv.org/abs/2507.23415) already proposes EB-induced
  spin magnetization, including intrinsic response in time-reversal-symmetric
  insulators. The mere existence of such a response is not our novelty.
- [Xiao et al., npj Computational Materials 8, 138 (2022)](https://www.nature.com/articles/s41524-022-00828-1)
  derives nonsynchronous polarization/BPVE switching in sliding bilayers from
  their different tensor transformations. It is close conceptual prior art
  for “polarization does not track every response,” but concerns optical BPVE,
  not this EB spin coefficient or the present local power comparison.
- [Zheng et al., Nano Letters (2025)](https://arxiv.org/abs/2402.17166) identifies
  layer coherence as an origin of planar Hall effects tunable by sliding/strain.
  It is mechanism background rather than the same spin observable; “interlayer
  coherence matters” alone is already established.

These sources constrain claims; they are not an exhaustive proof of priority.
Before promoting this to a paper direction, resolve the driven/static mapping,
derive transferable symmetry conditions, and verify a useful magnitude in a
realistic model/material. No DFT or material-specific claim is made here.
