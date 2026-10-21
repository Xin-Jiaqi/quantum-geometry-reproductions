# B+D Gate 5 | The interlayer spin-transfer / local-torque cancellation gate

**Date:** 2026-10-21. **Independent toy physics; not Xiang/Jin/Wang PRL (2025) paper reproduction.**  
**Status:** MODEL WARD CHECK PASS; observable material response / new mechanism **OPEN**. No DFT, finite-width edge, spin-lattice dynamics, relaxation/bath, or real-ferroelectric demonstration.

Executable: [32-band spin-transfer toy](../scripts/bd_gate5_torque_ward.py), which imports the [Gate 4 32-band model](../scripts/bd_gate4b_atomic_spin.py) and honeycomb geometry. \`python scripts/bd_gate5_torque_ward.py --dense\` runs meshes 11, 19, 23. The test targets a *different* physical question than the static in-plane spin-density tables in Gate 4.

## Idea

Could real polar magnetic stacking produce a **large layer-to-layer spin-transfer torque** that is hidden by a near cancellation with local exchange/SOC torque, leaving a much smaller net layer-resolved spin accumulation? This might motivate a frequency-resolved probe of the **internal angular-momentum pathway** rather than only a net Edelstein moment.

**This is a screening hypothesis, NOT an original mechanism claim**: spin continuity, torque partition and torque cancellation are standard. We would need a distinct *microscopic interlayer-resolved observable*, robustness and an experiment that cannot be accounted for by insulating Néel spin-orbit torque or ordinary spin magnetoelectricity.

## Operators and exact Ward identity

In the 32-band spin-independent spatial hopping + atomic SOC + A-type-AFM model, choose the physically local layer-staggered spin density

$$
S_-^x=\frac12 \tau_z\sigma_x,
\qquad H=H_{\mathrm{inter}}+H_{\mathrm{local}}.
$$

Define the two **operator-exact** torques, not an arbitrary current partition:

$$
\mathcal T_{\mathrm{inter}}^x=i[H_{\mathrm{inter}},S_-^x],
\qquad
\mathcal T_{\mathrm{local}}^x=i[H_{\mathrm{local}},S_-^x].
$$

Because \(\hbar=1\),

$$
\dot S_-^x=\mathcal T_{\mathrm{inter}}^x+\mathcal T_{\mathrm{local}}^x.
$$

The retarded length-gauge occupied–unoccupied Kubo result in the *zero-broadening subgap* independent-electron model obeys

$$
\boxed{\chi_{\mathcal T_{\mathrm{inter}},E_x}(\omega)+
\chi_{\mathcal T_{\mathrm{local}},E_x}(\omega)
=-i\omega\chi_{S_-^x,E_x}(\omega)}
$$

after all bands are included, with the tested gauge/contact prescription. This is an identity, **not a newly derived universal geometric coefficient**. Here \(H_{\mathrm{inter}}\) is defined unambiguously as \`model(out_scale=1) - model(out_scale=0)\` at fixed structural registry and all other parameters; the local part contains the intralayer spin-independent hopping, onsite SOC and exchange. In a material, partition-dependent spin fluxes and torque dipoles require further evaluation.

## Numerical results and actual convergence

Physical registry: AB fractional shift \((1/3,1/3)\) (see Gate 4), \(m=0.34\), atomic SOC \(=0.30\), synthetic \(\omega=0.06\), \(\eta=0\). The direct gap is \(\sim0.12356\), so this is **strictly below the interband threshold**. Units and normalizations are synthetic; the response is a *uniform electric-drive susceptibility* not a measured flux.

| Mesh | \(\chi_{S_-^x,E_x}\) real | \(\mathrm{Im}\chi_{\mathcal T_{\rm inter},E_x}\) | \(\mathrm{Im}\chi_{\mathcal T_{\rm local},E_x}\) | \(\mathrm{Im}\chi_{\mathcal T_{\rm total},E_x}\) |
| ---: | ---: | ---: | ---: | ---: |
| \(7^2\) | \(+3.668962\times10^{-6}\) | \(+8.707718\times10^{-5}\) | \(-8.729731\times10^{-5}\) | \(-2.201377\times10^{-7}\) |
| \(11^2\) | \(+8.422135\times10^{-6}\) | \(+9.728828\times10^{-5}\) | \(-9.779361\times10^{-5}\) | \(-5.053281\times10^{-7}\) |
| \(15^2\) | \(+9.632840\times10^{-6}\) | \(+1.017808\times10^{-4}\) | \(-1.023588\times10^{-4}\) | \(-5.779703\times10^{-7}\) |
| \(19^2\) | \(+9.696441\times10^{-6}\) | \(+1.021557\times10^{-4}\) | \(-1.027375\times10^{-4}\) | \(-5.817865\times10^{-7}\) |
| \(23^2\) | \(+9.698690\times10^{-6}\) | \(+1.021742\times10^{-4}\) | \(-1.027562\times10^{-4}\) | \(-5.819214\times10^{-7}\) |

At \(23^2\), the total torque response magnitude is only \(\simeq0.5695\%\) of the interlayer-transfer contribution. The Ward residual is \(\sim7.3\times10^{-18}\) at \(23^2\); analytic operator continuity residual is at floating-point noise. At \(\omega=0\), net torque vanishes as expected; at finite subgap \(\omega\) a nonzero total rate of change remains.

**Switching and disconnected controls:** with the same A-type AFM \(m\), the AB \(\leftrightarrow\) BA states related by \(m_zT\) have opposite layer-staggered density and both torque susceptibilities at *every tested mesh*. At \(11^2\), setting any of the interlayer hopping amplitude, atomic SOC, or magnetic exchange to **zero** suppresses the entire selected response to \(10^{-11}\) or below. For **total** \(S_{\mathrm{tot}}^x=\sigma_x/2\), the direct interlayer torque is exactly zero because the interlayer Slater–Koster hopping is spin independent; layer-staggered transfer is a distinct *operator-level bookkeeping channel*, not proof of spin transfer into the magnet.

The \(\omega=0.03,0.06,0.09\) values independently reproduce the Ward identity at \(11^2\), but **frequency-dependent broadening/scattering, finite q, magnetic resonance, torque-dipole and experimental readout are not modeled**. No magnon excitation energy is assigned.

## Scientific interpretation and no-go boundaries

1. **Mathematically real:** intralayer/local spin torque and layer-to-layer transfer can be much larger than net \(\dot S_-^x\) in this toy model. The near cancellation and connection sign have explicitly converged.
2. **Not proved as new physics:** this cancellation is broadly required by continuity when local spin density is small/slow; it does not define an independent “spin quantum geometry” tensor or prove quantization.
3. **Not yet a magnetic-torque readout:** exchange torque \(i[H_{\rm exchange},S]\), SOC torque \(i[H_{\rm SOC},S]\) and the spin-transfer contribution must be individually tied to a **magnetic-order equation of motion** and separated from electronic spin precession before calling this a Néel switching torque.
4. **Not material-feasible yet:** Gate 4's toy AB/BA are metastable in its 32-band total-energy model, with nonpolar AA lower; actual real magnetic sliding-FE bilayer and magnetic texture still need verification.
5. **Not a certified Wannier spin code:** source Stage 1 top-K \`wannier90_hr.dat\` is insufficient without same-gauge spin, layer/projector, position and electromagnetic vertices.

## Direct literature competition to resolve before material jobs

- [Tang, Zhang & Cheng, *Nat. Commun.* 16, 7790 (2025)](https://doi.org/10.1038/s41467-025-63171-1) on insulating antiferromagnetic Néel spin-orbit torque and AC magnetic dynamics.
- [Shi et al., *PRL* 96, 076604 (2006)](https://doi.org/10.1103/PhysRevLett.96.076604) on physically proper spin current vs conventional current and torque dipole.
- [Zhou et al., *PRL* 137, 076401 (2026)](https://doi.org/10.1103/qhvs-pfkq) on layer Edelstein response.
- [*A unified symmetry framework for spin–ferroelectric coupling in altermagnetic multiferroics*, *Nat. Commun.* (2026)](https://doi.org/10.1038/s41467-026-69635-2) on multiple FE switching states and spin mapping in MnPS3.
- [Zhu et al., *PRL* 135, 056801 (2025)](https://doi.org/10.1103/dmzg-ck2t) on magnetic stacking connection operators.

**Priority decision:** retain this observable as **Gate 5 MODEL PASS / PHYSICAL NOVELTY OPEN**. On a real material, first calculate a *validated local exchange torque response* and frequency-resolved susceptibility in the **same gauge** as the structural FE pair. Do not commit to expensive DFT scans based on the present cancellation alone.
