# Round 2 | PT-staggered spin-geometry gate (2026-10-21)

**Status: conditional model-level PASS; novelty NOT established.** This report is independent post-paper brainstorming, **not** a result or correction of Xiang–Jin–Wang, PRL 135, 146303 (2025). It does not update the published paper's reproduction status.

## Question

Can a globally PT-symmetric magnetic bilayer have a **vanishing total** spin-metric-like interband kernel and a **nonzero layer-staggered** kernel, with a nonzero measurable *metallic* response after additional symmetry constraints are removed?

Here "spin metric" denotes the **real interband spin-current/velocity spectral product**, not a positive-definite Fubini–Study metric. Positivity and physical-current interpretation cannot be assumed.

## Exactly PT-symmetric four-band test

Natural units with independent layer Pauli matrices tau and physical-spin matrices sigma:

$
H(\mathbf k)=\varepsilon(\mathbf k) I+t\tau_x+\alpha\tau_z(\sin k_y\sigma_x-\sin k_x\sigma_y)+m\tau_z\sigma_z,
$
$
\varepsilon(\mathbf k)=-2b(\cos k_x+\cos k_y)+\gamma\sin k_x.
$

`PT = (tau_x i sigma_y) K`, acting at fixed k, with `(PT)^2=-1`. The last scalar term is a **synthetic symmetry-allowed PT term**, *not* derived from a realistic sliding geometry. For a magnetic crystal it requires an explicit microscopic hopping/flux explanation before being used in a material claim. In this test `t=.65, alpha=.55, m=.72`.

The two twofold-degenerate manifolds have spectral projectors `P_-` and `P_+`; their gap is `Delta=2 sqrt(t^2+m^2+alpha^2(sin^2 k_x+sin^2 k_y))`. Do **not** use independently gauged Abelian formulas inside exact PT doublets.

Define the Hermitian conventional spin current `J={v_x,sigma_y}/4`, and **one chosen** staggered-current convention `J_stag={tau_z,J}/2`. The U(2)-invariant interdoublet kernels are

$
G_O(\mathbf k)=\frac{\operatorname{Re}\operatorname{Tr}[P_-(\mathbf k)O(\mathbf k)P_+(\mathbf k)v_x(\mathbf k)]}{\Delta(\mathbf k)^2}.
$

The antiunitary PT transformation makes `J` odd and `J_stag` even, while `v_x` and `P_±` are even. Accordingly, `G_J(k)=0` but `G_Jstag(k)` is symmetry-allowed. **This is standard representation theory; it is not a new symmetry theorem.**

## Independent numerical tests

Script: [pt_hidden_metric_and_drude.py](../scripts/pt_hidden_metric_and_drude.py). Python NumPy; run directly.

At k=(0.8,0.6), `G_Jstag=-0.0872046915191`; `G_J` is zero to floating-point noise. Explicit checks pass: `(PT)^2=-1`, `PT H PT^{-1}=H`, `P^2=P`, PT parity of both operators, and invariance under **independent random U(2) rotations of the two doublets**.

51×51 periodic Brillouin-zone average:

| gamma | total G | staggered G |
| ---: | ---: | ---: |
| 0 | ~0 | ~0, despite nonzero local G |
| +0.15 | ~0 | -0.008420775 |
| +0.31 | ~0 | -0.017402935 |
| -0.31 | ~0 | +0.017402935 |

For `gamma=0` an extra momentum cancellation removes the integrated signal. For nonzero gamma the PT-required total zero survives while the staggered kernel survives integration. The gamma reversal is an **internal synthetic-model comparison, not a demonstrated BSF switch**.

### Transport reality gate

A nonzero integrated interband `G` is **not itself** a dc spin current, especially in an insulator. Using the explicitly chosen operator and a constant-tau *diagonal occupation RTA*, evaluate the Fermi-surface-only numerical coefficient

$
C_O=\int_{\mathrm{BZ}}\!\frac{d^2k}{(2\pi)^2}\sum_{n\in\mathrm{doublets}}[-f'(\varepsilon_n)]v_{n,x}\operatorname{Tr}(P_n O).
$

This convention keeps the doublet degeneracy through the projector trace; dimensional electric-charge/tau factors have been suppressed. It is a model-only Drude-type check, **not the original paper's full linear spin-response susceptibility**.

At `mu=0,T=.15,b=.72` (metallic, negative indirect gap), `C_J=0` and `C_Jstag` is:

- gamma=0: approximately 0;
- gamma=+0.31: `-0.03048982`;
- gamma=-0.31: `+0.03048982`.

For an actual indirect-gap **insulator** `b=.12, gamma=.31, mu=0, T=.02` (gap `0.81179445`), the spectral `G_Jstag=-0.017402935` is still nonzero, while `C_Jstag=-5.11e-12` is numerically negligible. **Thus an integrated "metric" does not guarantee a DC insulator Drude signal.** Additional contact/operator-derivative terms in the literal original SM must not be silently dropped when extrapolating to lattice Hamiltonians.

## What this does *not* prove

1. A new measurable "hidden spin quantum metric" effect. PT-staggered/two-layer cancellation and hidden geometric responses have numerous precursors.
2. That `J_stag` is the unique physically meaningful layer-resolved spin-flow operator. Spin is generally not conserved with SOC; site/layer continuity and boundary torque must be handled.
3. That `gamma sin k_x I` occurs along a dynamically stable ferroelectric sliding path or has a microscopic hopping origin in the proposed stacking systems.
4. That the metal RTA result survives realistic vertex corrections, finite-q edges, or magnetization dynamics.
5. A claim of an insulating Drude effect. The insulator test specifically rules that out in this toy model.

## Priority/novelty boundary: directly relevant external articles

- Xiang, Jin, Wang, [PRL 135 146303 (2025)](https://doi.org/10.1103/14mp-263q): original spin QGT and extrinsic linear spin-current geometry.
- [Sliding-ferroelectric layer spin Hall (J. Phys. Chem. Lett. 2024)](https://doi.org/10.1021/acs.jpclett.4c01479): sliding, layer-locked spin Berry curvature and spin Hall.
- [Layer spin Hall in altermagnets (PRB 112 134421, 2025)](https://doi.org/10.1103/643d-wkc1): layer spin Hall already considered.
- [Quantum geometric engineering in antiferromagnetic bilayers (Adv. Sci. 2025)](https://doi.org/10.1002/advs.202505860): magnetic coupling changes Hall/metric mechanisms.
- [PT sliding quantum-metric nonlinear valley Hall (PRB 113 235116, 2026)](https://doi.org/10.1103/p9sb-k9fj): PT-preserving hidden quantum geometry under sliding already demonstrated.
- [Stacking theory for 2D magnets (PRB 113 214455, 2026)](https://doi.org/10.1103/g471-xrkl): general magnetic stacking theory overlaps any proposed purely structural classification.
- [Sliding ferroic tuning of spin Hall in WTe2 (arXiv:2606.22015)](https://arxiv.org/abs/2606.22015): precise recent competition for an FE-switchable spin Hall claim.
- [PT-AFM hidden BCD nonlinear layer Hall (arXiv:2510.23971; published PRB 2026)](https://arxiv.org/abs/2510.23971): demonstrates the importance of distinguishing an integrated hidden geometry from experimentally detected nonlinear signal.

These are **abstract-level / accessible-text screening**, not complete novelty clearance. No claim of priority is made.

## Next decisive gates

1. Construct `H(k,delta)` from explicit geometry-dependent **microscopic interlayer hoppings** and real BSF connection operators. The two states must be distinct *physical structures*, not basis rephasings or sign changes in unrelated magnetic/hopping parameters.
2. Show how **layer-resolved spin density, bulk spin flow, boundary flux, and local SOC/exchange torque** close one continuity equation. If the "hidden" metric contributes only to a convention-dependent spin-current partition, discard as an independent discovery.
3. Use non-Abelian doublet projectors throughout and compare **full response** vs the isolated interband kernel, at both zero and finite omega.
4. Check existing layer-spin-Hall / hidden-geometry studies at equation/operator level before applying DFT to a candidate material.

**Interim verdict:** mathematically consistent, physically interesting as a screening diagnostic; novelty and device observability **OPEN**.
