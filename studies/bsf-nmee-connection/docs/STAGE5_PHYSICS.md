# Stage 5: weak-coupling coefficient, parity competition and fixed-stacking sign reversal

## What is actually established

We study an **eight-band synthetic bilayer Kane–Mele model**, not an actual ferroelectric sample. The observable is the **zero-temperature, gapped equilibrium static** coefficient obtained by differentiating the occupied-to-empty mixed spin–momentum Berry curvature with respect to a weak perpendicular Zeeman field. It has *not* been identified with the unresolved driven-DC NMEE definition of Jia et al., PRL 136, 226302 (2026).

The more precise physical mechanism is **competition**, rather than a claim that the interference term reverses sign during sliding. For the frozen comparison below, the response is negative at small sliding because a negative *odd-in-relative-cross-hopping-sign* contribution dominates. At larger sliding, a positive *even* contribution grows and overcomes that negative contribution. Both pieces are defined as controlled comparisons of fully recomputed Hamiltonians in an explicit orbital basis; only their sum is observable.

## Two hopping channels and symmetry-enforced weak-coupling expansion

At any fixed sliding vector, decompose the interlayer hopping into two matrices in a **specified, fixed orbital basis**:

\[
H(k; s,c)=H_0(k)+s V_{\mathrm{same}}(k)+c V_{\mathrm{cross}}(k).
\]

Each `V` contains both upper-right and lower-left hopping blocks and their correct phase/velocity derivatives. Here `same` connects A-to-A and B-to-B, while `cross` connects A-to-B and B-to-A. These labels are basis dependent. In the actual code, `s` and `c` are independent hopping energies, implemented by `same_scale`/`cross_scale` times `tunneling_eV`.

The layer-sign unitary `W=diag(I_4,-I_4)` maps `H(s,c)` to `H(-s,-c)` and leaves the static total spin-electric response unchanged. If an insulating gap stays open and the response is analytic, the Taylor expansion must therefore have **even total degree**:

\[
\Gamma(s,c;u)=\Gamma_0(u)+ A(u)s^2+B(u)c^2+C(u)sc+O((|s|+|c|)^4).
\]

For this model, `Gamma_0=0` because the two decoupled layers have mirror-opposite Rashba spin–momentum responses. The leading mixed (interference) coefficient is uniquely defined **within the fixed hopping-channel parameterization** by

\[
C(u)=\left.\frac{\partial^2\Gamma}{\partial s\partial c}\right|_{s=c=0}
=\lim_{t\to0}\frac{\Gamma(t,t;u)-\Gamma(t,-t;u)}{2t^2}.
\]

The second equality relies on the exact layer-sign equivalence. The completely independent four-quadrant central difference agrees. This is a genuine well-defined mixed *perturbation derivative*, but not an invariant statement about splitting a particular band-resolved geometric quantity without specifying the orbital gauge and physical perturbations.

For the occupied spectral projector `P(k,h_x,B_z;s,c)` of `H-h_x S_x-g mu_B B_z S_z`, the underlying static mixed Berry curvature is

\[
\Omega_{k_x,h_x}=i\operatorname{Tr}\big(P[\partial_{k_x}P,\partial_{h_x}P]\big),
\qquad
\Gamma\propto\left.\partial_{B_z}\int_{\mathrm{BZ}}\Omega_{k_x,h_x}\,\frac{d^2 k}{(2\pi)^2}\right|_{B_z=0}.
\]

Thus the leading mixed coefficient can be expressed **basis-independently for two specified physical perturbations** as

\[
C=\left.\partial_s\partial_c\partial_{B_z}\int_{\mathrm{BZ}}i\operatorname{Tr}\big(P[\partial_{k_x}P,\partial_{h_x}P]\big)\,\frac{d^2 k}{(2\pi)^2}\right|_{s=c=B_z=0}.
\]

To see the virtual-hopping origin, use the Riesz occupied projector `P=(1/(2 pi i))∮ (z-H)^(-1) dz`, where `G=(z-H0)^(-1)`:

\[
P_s=\frac{1}{2\pi i}\oint GV_{\mathrm{same}}G\,dz,\qquad
P_{sc}=\frac{1}{2\pi i}\oint (GV_{\mathrm{same}}GV_{\mathrm{cross}}G+GV_{\mathrm{cross}}GV_{\mathrm{same}}G)\,dz.
\]

This is a rigorous formal representation of the **two virtual interlayer hopping insertions**, with all occupied/empty subspace mixing included. It does not, by itself, fix the sign or imply universality. A symmetry selection rule for nonzero `C` must transform BOTH perturbation matrices and the response tensor; we have not proven such a rule for a real layer group.

## Numerical tests: the leading coefficient is negative at all three checked slides

The values below come from `results/sliding_stage5_weak_coupling.json`, with a parabolic fit of `Gamma/t²` against `t²` at `t=1.25,2.5,5 meV`, and valley-patch integration. Units of A, B, C are `response units / eV²`, **not experimental response units**.

| Sliding u | leading total A+B+C | leading even A+B | leading mixed C |
|--:|--:|--:|--:|
| 0.10 | -0.02680 | +0.003148 | **-0.02995** |
| 0.20 | -0.03284 | +0.023978 | **-0.05682** |
| 0.40 | +0.07677 | +0.168239 | **-0.09147** |

**Correction to stage-4 interpretation:** the stage-4 *inclusion-exclusion nonadditivity* at `t=20 meV` changes sign with slide, but that does **not** mean the leading `sc` interference coefficient changes sign. Indeed `C<0` at all three checked shifts. Inclusion-exclusion at finite t includes higher-order contributions such as `s² c²`, which are even under `c → -c`. In contrast, the odd-in-c parity component isolates terms odd in both s and c, including `sc, s³c, sc³, ...`.

At `t=20 meV` (numbers in code response units):

| Slide u | even in cross-sign | odd in cross-sign | total |
|--:|--:|--:|--:|
| 0.10 | +6.3379×10^-6 | -23.3187×10^-6 | **-16.9808×10^-6** |
| 0.20 | +38.5041×10^-6 | -37.5719×10^-6 | **+0.9322×10^-6** |
| 0.40 | +138.3592×10^-6 | -39.8571×10^-6 | **+98.5021×10^-6** |

Therefore the model sign reversal is a **competition between a positive even channel and negative relative-phase-sensitive interference**, rather than an interference coefficient switching from negative to positive. Note that `even/odd` is defined by comparing distinct Hamiltonians obtained by reversing the sign of *all* cross-sublattice interlayer hoppings while holding monolayers and same-sublattice hoppings fixed. This is not an experimentally demonstrated ferroic operation or an invariant decomposition under arbitrary orbital-basis changes.

## A distinct result: fixed stacking changes sign with *tunneling*, without gap closure

At **the same** stacking `u=0.20` and fixed SOC and chemical potential, increasing the common scale `t_perp` produces:

| `t_perp` | Gamma_xz;x | min sampled gap from mu |
|--:|--:|--:|
| 10 meV | -2.8303×10^-6 | >9.8 meV |
| 15 meV | -4.3567×10^-6 | >6.5 meV |
| 20 meV | +0.9322×10^-6 | >3.3 meV |

Valley-patch bisection (`61×61` exterior midpoint points, polar `60×80` per valley) finds a **fine-mesh root near `t_perp ≈ 19.63385 meV`**; the gap to `mu` remains ~3.59 meV. An additional 46×46 reciprocal-grid search with 14 local optimizations per band edge gives an **indirect valence-to-conduction gap ≈17.04 meV** at the root (and ≈16.76 meV at 20 meV tunneling). These are extensive numerical checks, not a mathematical proof of a global spectral gap. Repeating B derivatives with `0.125 T` rather than `0.25 T` and using the finer k mesh shifts Gamma by order `10^-9` or less near the root. The exact root depends on the model parameters and remains subject to numerical quadrature; this is **not** an equilibrium structural phase transition.

The fixed-slide sign reversal demonstrates that even holding the sliding displacement constant, higher-than-quadratic hopping processes can overcome the weak-coupling sign. The relative interference contribution remains negative in the checked range.

## Numerical limitations

- The coarse exterior midpoint BZ quadrature does not integrate a constant exactly when circular K/K' disks are excised. At the 41×41 exterior grid, **area-weight deficit is about 0.028 out of 19.088** (0.15%). At 61×61 it falls to ~0.001. Valley-concentrated response convergence tests yield much smaller changes in Gamma; this is a numerical finding, not a formal bound on arbitrary smooth integrands.
- Uniform large k meshes previously produced misleading fluctuations; polar valley patches and an explicit coarse-to-fine response check are mandatory.
- “Insulating” means sampled H(k) retains the chosen chemical potential within the gap, not a globally rigorous certificate at every continuum k.
- SOC, lattice embedding, orbital phases, the hopping range, and the spin/position operators remain synthetic choices. The model has no first-principles equilibrium stacking energies or ferroelectric polarization.
- The PRL intrinsic *driven-DC* vs full *static equilibrium* definition remains open pending correspondence with the authors. Every reported number here uses the **static** branch only.

## Connecting to a real material: one step at a time

We do **not** jump straight to high-throughput DFT. The immediate “material bridge” is to ask how `H0`, `V_same`, and `V_cross` would be obtained from a single **consistent Wannier orbital manifold**:

1. Choose a simple nonmagnetic, insulating sliding bilayer for which two symmetry-connected candidate registries and an interlayer spacing can be specified; the group toolkit gives their symmetry and connection operator, **not** their energy minima.
2. Freeze SOC-including monolayer orbital character and use a **common set of localized orbital projectors** for both layer blocks. A chosen Wannier function label or `A/B` sublattice is not universally meaningful: in real materials, grouping tunneling into symmetry/distinct orbital channels must be declared explicitly.
3. Obtain layer-block `H11(k)`, `H22(k)` and `T(k; tau)` from spinor Wannier models at **each relaxed registry**; compare symmetry-covariant, aligned orbital gauges before subtracting hopping terms. Slipped structures generally cannot be modeled simply by reusing a fit at one registry.
4. Supply full spin matrix and position/Berry-connection information consistently with the Wannier gauge. `wannier90_hr.dat` alone does **not** establish a quantitative NMEE. Check spectral gaps, k-mesh and Zeeman derivative and both response protocols.
5. Only after (1)-(4), compare decoupled reference, coupled response, and symmetry-partner response. No real-material prediction has been obtained in this stage.

## Decision about scientific novelty

**Positive:** the model now has a controlled weak-coupling mixed derivative, a clearly specified competition mechanism, a sign reversal at fixed sliding when tuning tunneling, and independent projector/finite-Zeeman checks. The deeper question is whether orbital symmetries *force* or *suppress* the `sc` coefficient across a class of BSF stackings. **Not yet proven:** a material-independent sign rule, a real switchable ferroelectric bilayer with sizable total response, or full equality with the PRL driven DC NMEE. PRL-level claims remain premature.
