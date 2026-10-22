# BD25 test 1 | What the published LEE minimal model can and cannot test

**Date:** 2026-10-22  
**Status:** BASELINE_NUMERIC_TEST_PASS / TARGET_RESPONSE_NO_GO_WITHIN_EQ3 / MATERIAL_NOVELTY_UNRESOLVED  
**Compute:** zero DFT, zero HPC.

## 1. The honest result

The proposed target was the current-induced layer-odd out-of-plane spin response
\[
\chi_{S_-^z,j_y}(u_x),\quad S_-^z=S_T^z-S_B^z.
\]
Before adding any sliding term, I audited the source Hamiltonian in Zhou et al., [*Layer Edelstein Effect*, PRL 137, 076401 (2026), Eq. (3)](https://arxiv.org/html/2607.17209):
\[
H(\mathbf k)=\epsilon_0+t_\perp\tau_x+\Delta\tau_z+
\alpha\tau_z(\sigma_xk_y-\sigma_yk_x).
\]

Define the in-plane spin-helicity operator
\[
\Sigma_{\mathbf k}=(k_y\sigma_x-k_x\sigma_y)/k.
\]
It commutes with this Hamiltonian. In a helicity sector \(s=\pm1\), the spin state is an eigenstate of an **in-plane** spin operator and the remaining \(2\times2\) Hamiltonian acts only in layer space:
\[
H_s=\epsilon_0+t_\perp\tau_x+(\Delta+s\alpha k)\tau_z.
\]
Thus each nondegenerate band can be chosen as a product of a layer state and an in-plane spin-helicity state, with \(\langle\sigma_z\rangle=0\). Consequently
\[
\boxed{\langle\tau_z\sigma_z\rangle_{n\mathbf k}=0}
\]
for all \(k,\Delta,t_\perp,\alpha\). In the paper's **band-diagonal Boltzmann/Eq. (2) regime**, the target response is therefore identically zero in Eq. (3). Replacing \(t_\perp\) or \(\Delta\) by scalar functions of registry \(u\) does not change this result.

**This is a useful model-level NO-GO, not a no-go theorem for real bilayers.** The target invariant allowed by our earlier \(C_{3z}+m_x+m_z\) symmetry calculation cannot be numerically tested using Eq. (3) alone. A source-backed \(S_z\)-sensitive SOC/valley/interlayer matrix element is missing. I will not invent a spinful hopping just to force a nonzero \(\gamma\).

## 2. Baseline response curve from the published model

A NumPy-only script evaluates the same radial Boltzmann structure as the paper's Eq. (2), for the existing in-plane responses. This is a dimensionless **sanity calculation**, not a fit to MoSSe/MoTe₂/WTe₂ and not a reproduction of the DFT curves.

Parameters: \(t_\perp=1,\alpha=0.5,m=1,\mu=0,T=0.02\), \(k_{\max}=4.5\), 16,000 radial points. The response is normalized by the charge-current weight; common units/constant factors are suppressed.

| \(\Delta/t_\perp\) | Layer-odd \(\chi_{S_-^y,j_x}\) | Global \(\chi_{S_+^y,j_x}\) | Layer-odd \(\chi_{S_-^z,j_y}\) |
|---:|---:|---:|---:|
| -0.8 | 0.31070 | +0.11634 | 0 |
| -0.4 | 0.42452 | +0.08182 | 0 |
| 0 | 0.48497 | 0 | 0 |
| +0.4 | 0.42452 | -0.08182 | 0 |
| +0.8 | 0.31070 | -0.11634 | 0 |

Parity tests on 161 points passed: layer-odd in-plane response is even in \(\Delta\); global in-plane response is odd in \(\Delta\); \(S_z\) is exactly zero in this minimal model. The \(S_y\) trends are baseline sanity checks consistent with the model's layer compensation/gate unbalancing, not new physics.

Code: [\`lee_published_model_baseline.py\`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/lee_published_model_baseline.py). It writes a CSV and repeats the parity assertions when run.

## 3. What this teaches us about the idea

1. The earlier group-theory result remains correct: symmetry allows one first-order invariant \((u_xj_y-u_yj_x)S_-^z\) for the **idealized** group and specified layer-exchange representation.
2. The allowed coefficient is not fixed by symmetry and is zero in the published Eq. (3) model because that model has no out-of-plane spin texture.
3. The full MoSSe valley model described in the user's reading notes contains spin-mixing terms, but the source does **not** provide a registry-dependent interlayer Hamiltonian/force vertex. We still cannot infer \(\partial_{u_x}\chi\) from it.
4. The 2026 LEE paper is the **closest framework competitor**, but the accessible main-text model treats \(E_z\)-controlled layer detuning and discrete stackings, not this continuous slide derivative. The supplemental equations must be compared before deciding whether the exact derivative is already included. Do not declare either novelty or pre-emption yet.

## 4. Next physical gate

Before any \(\chi(u)\) curve, obtain one **source-backed** model with both (i) a nonzero \(S_z\) spin texture or spin-z matrix elements and (ii) a microscopic registry derivative \(\partial_u H\), such as a published spinful bilayer tight-binding/continuum model or source-verified Wannier Hamiltonian. Then:

- derive \(\partial_u H\) from its stated hopping/orbital geometry rather than choosing an arbitrary coefficient;
- calculate \(\chi_{S_-^z,j_y}(u_x)\) and test \(\chi(u_x)=-\chi(-u_x)\);
- include a mirror-restored \(u=0\) control and an SOC-off control;
- compare the exact observable, drive and transport protocol with the 2026 LEE paper before assigning literature role.

If no such source model can be found cheaply, the rational move is to retire **this specific S_z slide seed as a computable candidate**, preserve its correct symmetry result, and choose a different response already present in a source-backed Hamiltonian. No DFT/Slurm/WSL/Japan task is authorized.


## 5. A viable source-backed model route found after the no-go

A plausible next model is **aligned bilayer MoS\(_2\) in the high-symmetry AA registry**, not a magnetic domain wall. The monolayer 1H structure has a horizontal mirror; the AA bilayer is the natural candidate for a layer-swapping \(m_z\) reference. This must be confirmed on the actual atomic coordinates before use.

- Venkateswarlu, Honecker & Trambly de Laissardière, [*Electronic localization in twisted bilayer MoS\(_2\) with small rotation angle*, PRB 102, 081103(R) (2020)](https://perso.u-cergy.fr/~ahonecker/papers/tmos2-tbPRBRfinal.pdf), supplies an 11-orbital-per-layer Slater–Koster tight-binding model and **registry-sensitive interlayer hopping** including S–S, Mo–S and Mo–Mo terms. It was benchmarked against DFT near the gap, especially valence bands; the paper warns that its relaxed-structure effects are not fully included and its small-angle results are qualitative.
- A spinful extension must use a separately source-backed atomic SOC parameterization, e.g. Zhan et al., [*Multi-ultraflatbands tunability and effect of spin-orbit coupling in twisted bilayer transition metal dichalcogenides*, arXiv:2005.13868](https://arxiv.org/abs/2005.13868). This is a model-construction source, not yet a validated aligned-AA implementation.
- The recent [*Magnetization-dependent and stacking-tunable Edelstein effect in the two-dimensional magnet*, PRB 114, 125129 (2026)](https://doi.org/10.1103/48gx-cpg2) and the 2024 magnetic-sliding-ferroelectric CrI\(_3\) work show that the **broad claim** “sliding can tune a current-induced Edelstein response” is not new. Their magnetic observables/materials are not identical to the proposed nonmagnetic layer-odd \(S_z\) derivative, so they are **Introduction/strong-neighbor context** for that narrower claim, not automatic direct pre-emption.
- The 2020 SK model alone is not enough: its published implementation does not establish the desired spin response. Before any full \(k\)-mesh, first verify (i) the AA coordinate-level \(m_z\), (ii) the spinful SOC Hamiltonian, (iii) the finite-difference derivative of the **actual** interlayer hopping under \(u\), and (iv) that the model gives the correct spin-split monolayer/bilayer bands.

**Proposed next action:** build only the real-space hopping/structure unit test for two unrelaxed aligned layers at \(u=0,\pm\delta\), and confirm the layer-swap covariance \(H(u)=U_{m_z}H(-u)U_{m_z}^{-1}\). No spin-response integration until this passes. This is a source-backed path to \(\partial_uH\), unlike adding a free \(\tau_z\sigma_z(\mathbf u\times\mathbf k)_z\) term. It remains a moderate implementation task and has **not** been coded or run yet.


### 5.1 Geometry gate now implemented and run

The next tiny gate has been completed in [\`aa_mos2_registry_geometry_gate.py\`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/aa_mos2_registry_geometry_gate.py), using the source geometry above: two unrelaxed aligned AA layers, Mo–Mo separation \(6.80\) Å, lattice constant \(3.18\) Å, S height \(\pm0.49115a\), and the explicit internal S-upper/S-lower swap under \(m_z\).

**Local test result:** 450 atom-pair/lattice-translation comparisons for two nonzero displacements passed \(d_{ab}(\mathbf u)=d_{\bar b,\bar a}(-\mathbf u)\), with maximum absolute distance mismatch \(1.776\times10^{-15}\) Å. The interfacial S–S plane spacing is \(3.6763\) Å, consistent with the source model's tabulated \(3.676\) Å. The initial test exposed a code bug in the species-label assertion (the two chalcogen planes are the same chemical species); this was corrected before the passing run.

This is only the **geometric covariance prerequisite**. It does not prove the electronic hopping matrix covariance, because orbital parity/rotation matrices and the Slater–Koster angular factors have not yet been implemented. Next step: implement one spin-independent SK hopping block and verify the full matrix relation \(H(\mathbf u)=U_{m_z}H(-\mathbf u)U_{m_z}^{-1}\) before adding SOC or calculating response.
