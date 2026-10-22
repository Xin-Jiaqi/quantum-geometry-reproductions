# B+D Gate 6 | What actually cancels the interlayer electronic spin-transfer torque?

**Study date:** 2026-10-22. **Status: independently reproduced SYNTHETIC MODEL DECOMPOSITION, NO NOVELTY CLAIM.** No HPC, new MnPS3 VASP, R2 submission or material properties are involved. Original Xiang/Jin/Wang (2025) PRL reproduction and its acceptance are untouched. This is a **re-audit of the same Gate 4–5 32-band toy Hamiltonian**, not a new fitted Hamiltonian.

**Executable:** [\`bd_gate6_exchange_soc_torque.py\`](../scripts/bd_gate6_exchange_soc_torque.py), importing the existing Gate 4A and Gate 4B model modules. An independent local reconstruction of the same microscopic 32-band model gave the dense-mesh values below before archival; the committed script recomputes those values with an independent decomposition directly from the existing Gate 4 model.

## 1. Scientific question and operator choice

Gate 5 found an apparent giant layer-transfer torque nearly compensated by local torque. **Which local interaction cancels it?** If this is mostly atomic SOC, it is not necessarily a torque on the Mn local moments. If it is magnetic exchange, a magnetic angular-momentum sink is possible — but electronic-exchange torque **still does not by itself compute Néel switching**.

For the **existing A-type interlayer-AFM toy**, use the same layer-staggered electronic spin density \(S_-^x=\tau_z\sigma_x/2\). The Hamiltonian already has exactly:

- \(H_\perp\): spin-independent interlayer Slater–Koster hopping;
- \(H_{\mathrm{intra,hop}}\): spin-independent intralayer hopping;
- \(H_{\mathrm{ex}}=m\tau_z\sigma_z\): fixed classical layer-staggered exchange, \(m=0.34\);
- \(H_{\mathrm{SOC}}=\lambda_{\mathrm{SOC}}\mathbf L\cdot\boldsymbol\sigma/2\), \(\lambda_{\mathrm{SOC}}=0.30\).

With \(\hbar=1\), define \(\mathcal T_a=i[H_a,S_-^x]\). Since \(H_{\mathrm{intra,hop}}\) is spin independent and layer diagonal, its contribution is exactly zero; the **operator** identity is

$$
\boxed{\dot S_-^x=\mathcal T_\perp^x+
                 \mathcal T_{\mathrm{ex}}^x+
                 \mathcal T_{\mathrm{SOC}}^x.}
$$

In this toy representation specifically,

$$
\mathcal T_{\mathrm{ex}}^x=
i[m\tau_z\sigma_z,\tau_z\sigma_x/2]
=-m\sigma_y=-2mS_{\mathrm{tot}}^y.
$$

This exact algebra is a useful **operator identification**, not a new torque mechanism.

Crucially, the generalized force for a *hypothetical rigid rotation of the toy Néel axis about \(y\)* is a **different object**. If \(H_\mathrm{ex}(\theta)=m\tau_z(\sigma_z\cos\theta+\sigma_x\sin\theta)\), then

$$
F_\theta=-\frac{\partial H_\mathrm{ex}}{\partial\theta}\bigg|_0
=-m\tau_z\sigma_x=-2mS_-^x.
$$

For this simple constrained-spin parametrization, \(\chi_{F_\theta,E}=-2m\chi_{S_-^x,E}\). This **electronic generalized force** is a prospective input to an independently specified local-moment/magnon equation, not a magnetization dynamics calculation or automatically the full Néel torque.

## 2. Same full-band causal Kubo, now with four torque channels

Compute the **full occupied-to-empty interband** response of the same 32-band model at eight-electron-band filling, from exactly the existing velocity \(\partial_{k_x}H\), chemical potential in the model gap and \(\eta=0\), and \(\omega=0.06\) below minimum direct gap \(\approx0.1235575\). No Drude or fitted lifetimes. Use the original toy AB/BA pair, whose interlayer structural map is **\(m_zT\) in this toy**; this is not the actual MnPS3 \(L_2\) Eq.13.

| \(N_k\) | \(\Im\chi_{\perp,E_x}\) | \(\Im\chi_{{\rm ex},E_x}\) | \(\Im\chi_{{\rm SOC},E_x}\) | \(\Im\chi_{\rm net,E_x}\) | \(\Re\chi_{S_-^x,E_x}\) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| \(15^2\) | \(+1.0178080383\times10^{-4}\) | \(-1.0047647087\times10^{-4}\) | \(-1.8823033559\times10^{-6}\) | \(-5.7797039409\times10^{-7}\) | \(+9.6328399015\times10^{-6}\) |
| \(19^2\) | \(+1.0215574595\times10^{-4}\) | \(-1.0081096204\times10^{-4}\) | \(-1.9265703600\times10^{-6}\) | \(-5.8178645188\times10^{-7}\) | \(+9.6964408648\times10^{-6}\) |
| \(23^2\) | \(+1.0217424295\times10^{-4}\) | \(-1.0082850585\times10^{-4}\) | \(-1.9276584699\times10^{-6}\) | \(-5.8192137015\times10^{-7}\) | \(+9.6986895025\times10^{-6}\) |

All are dimensionless model coefficients. At \(23^2\):

- Exchange electron torque is **98.68%** of the interlayer term in opposite sign.
- SOC torque is only **1.89%** of that interlayer magnitude (also opposite sign).
- The **net** spin-rate signal is \(0.57\%\) of the interlayer term.
- The **generalized Néel-axis force** has \(\chi_{F_\theta,E_x}=-6.59510886\times10^{-6}\) (real at this subgap point), about **6.45%** of the interlayer torque response **magnitude**, but a different phase/physical observable.

The numerical electronic Ward residual is \(7.3\times10^{-18}\), and the separately assembled operator commutator is zero within machine noise. Relative results are nearly unchanged from \(19^2\) to \(23^2\). At \(11^2\) the toy AB/BA partner gives the sign-reversed layer-staggered responses and the separate component null controls \(t_\perp=0,\lambda_{\mathrm SOC}=0,m=0\) are zero within numerical precision. Frequency points \(\omega=0.02,0.06,0.09\) preserve the complete Ward relation, but **not** a low-\(\omega\) material resonance claim.

**Interpretation:** the "local torque" in Gate 5 is nearly entirely *electronic exchange torque* for this toy and only a small SOC contribution. The numerical cancellation is compatible with spin continuity but **the Ward identity does not, by itself, require the ratio 98.68%**. Neither that ratio nor the sign switch is a claim of an original quantum-geometric response.

## 3. Mandatory material-magnetism correction: layer-odd is NOT MnPS3 Néel-odd

The toy's \(m\tau_z\sigma_z\) describes an **A-type layer AFM**, with every site in a layer carrying the same local exchange sign. Real MnPS3 has a honeycomb **intrinsic Néel pattern within each Mn layer**. Its magnetic states \(L_1,L_2\) consist of **four Mn sites per bilayer magnetic cell** (two opposite-spin Mn sublattices in each layer; interlayer alignment further varies).

Let \(P_{\ell a}\) be **physical same-gauge Wannier/PAW site projectors** on layer \(\ell=T,B\), Mn sublattice \(a=A,B\), and \(\eta_{\ell a}=\pm1\) the **source-verified magnetic sign**. Distinguish

$$
S^i_{\mathrm{layer-odd}}=\sum_{\ell,a}l_\ell P_{\ell a}s^i P_{\ell a}
\quad(l_T=+1,l_B=-1),
$$

from a magnetically weighted Mn spin density

$$
S^i_{\mathrm{Mn\;N\acute eel}}=\sum_{\ell,a}\eta_{\ell a}P_{\ell a}s^iP_{\ell a}.
$$

These are *not the same operator*. For illustrative \(\eta=(+,-,+,-)\) or \((+,-,-,+)\) in site order (T,A;T,B;B,A;B,B), the simple layer weights \(l=(+,+,-,-)\) are orthogonal to \(\eta\). **These are illustrative magnetic patterns, NOT a claim about the actual verified source L₂ spin/site indexing.** The actual \(\eta\) must be imported from the private R1 structure preflight with its source site map.

The proper physical exchange torque on localized Mn moments is defined per site by a variation with **each local moment direction**, not by assigning the entire electronic layer torque to a single classical magnetic order parameter. An electronic spin-density/torque response may be sizable yet couple weakly to the *relevant Néel normal mode* after projection.

Therefore **the Gate 5/6 large exchange/interlayer cancellation cannot be transferred as a numerical or even qualitative MnPS3 prediction** without a new source-faithful sublattice exchange Hamiltonian and properly defined spin/position/site vertices.

## 4. Literature audit (fresh sources checked 2026-10-22)

- [Guo et al., *Layer-dependent spin-orbit torque switching of Néel vector in a van der Waals antiferromagnet*, *Nature Communications* **16**, 8911 (2025)](https://doi.org/10.1038/s41467-025-63966-2): real layer-dependent SOT switching in vdW A-type FCGT, a direct barrier to claiming "layer torque drives Néel order" is a new effect. **The material/magnetic texture differs from MnPS3**.
- [Haque et al., *Terahertz field-induced giant symmetry modulations in a van der Waals antiferromagnet*, *Nature Communications* (2026), online 20 October](https://doi.org/10.1038/s41467-026-78519-4): time-resolved SHG of **MnPS3**, phonon modes at **1.7 and 4.5 THz**, dynamical mirror-symmetry reduction. These are *driven lattice/phonon mechanisms*, **not** evidence for our undriven electronic spin-transfer/torque susceptibility. However plain "THz MnPS3" is pre-empted; separating electronic vs ionic channels is now mandatory.
- [Ilyas et al., *Coherent terahertz control of metastable magnetization in FePS3*, *Nature Communications* (2026), online 6 October](https://doi.org/10.1038/s41467-026-78164-x): competing THz spin–lattice driving in a related metal phosphorus trisulfide.
- [*Octupole-driven spin torque switching of antiferromagnetic tunnel junctions*, *Nature*, 21 October 2026](https://doi.org/10.1038/s41586-026-11172-5): clear precedent for relating inter-/intra-sublattice spin-current imbalance to AFM torque, though in metal AFM tunnel junctions (PtMn3), not a polar insulating MnPS3 bilayer.
- [Sun et al., *A unified symmetry framework for spin–ferroelectric coupling*, *Nature Communications* **17**, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2) and [Yang et al., *Symmetry-Programmed Pure Spin Photocurrents*, *ACS Nano* (2026)](https://doi.org/10.1021/acsnano.6c06799) already address MnPS3 FE magnetic-connection spin band / spin light-current switching. Our **full local-spin / exchange / ionic / spectral readout** question is distinct in definition but **novelty OPEN until their supplements and relevant SOT literature are matched equation-by-equation**.

## 5. Prioritized go/no-go physics gate (without another new toy)

**Gate 6a (PASS, model-scoped):** exact electronic transfer–exchange–SOC torque decomposition and a bounded eigenstate Kubo benchmark; independent material-general semantics corrected.

**Gate 6b (TO DO, read-only with private source structure):** reconstruct the published **L₂ magnetic four-Mn pattern** and its site projector representation; check actual SOC **magnetic** subgroup and classify the physical \(F_{\theta_\nu}=-\partial H/\partial\theta_\nu\) for each relevant Néel normal mode. Only then select \(\chi_{F_{\theta_\nu},E}(\omega)\); **do not assume \(S_-^x\) is a valid MnPS3 Néel mode**.

**Gate 6c (TO DO, independently validated solver):** in a complete site-resolved spin/position/Wannier or DFT matrix-element basis, calculate (i) full electronic response (ii) actual projected exchange force, (iii) an appropriate phonon/magnon contribution and experimentally meaningful detection metric, with correct boundary/vertex contacts and lifetime treatment. The real material R1 lacked even a converged ground-state SCF; NO physical numerical result exists.

**Decision:** Retain the idea "**ferroelectric registry changes the *partition and mode projection* of electronic angular momentum into magnetic and lattice channels**" only as a *falsifiable candidate*, not an already new mechanism. If it collapses to magnetic symmetry parity, conventional insulating NSOT or phonon-mediated THz control, archive as NO-GO before expensive production DFT.

**Operational:** \`tasks/CURRENT.md\` stays **R1_CLOSED; no auto jobs**. Any numerical SCF restart requires new R2 ID and user authorization. The pending next action is a **paper/site symmetry + source-derived mode projector review**, not submit VASP.
