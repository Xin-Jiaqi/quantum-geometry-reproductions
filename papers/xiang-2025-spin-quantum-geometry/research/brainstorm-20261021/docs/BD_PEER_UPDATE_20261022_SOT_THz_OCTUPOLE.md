# Four critical 2025–2026 competitors for B+D magnetic stacking physics

**Literature verified:** 2026-10-22 against the publisher's article pages.  
**Scope:** theoretical/experimental novelty screen for a *proposed* MnPS3 bilayer layer/site-resolved electronic spin-density, exchange-torque and finite-frequency response. No claim of a new material result; DFT R1 was closed without SCF convergence.  
**Main comparator for the candidate:** [Sun et al., *A unified symmetry framework for spin–ferroelectric coupling in altermagnetic multiferroics*, Nat. Commun. 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2).

## 1. Direct published papers (status, correct contrast and required falsification)

| Paper | Published and supported claim | Direct overlap with B+D | Distinction that must be *demonstrated*, not asserted |
| --- | --- | --- | --- |
| **Guo et al.**, [*Layer-dependent spin-orbit torque switching of Néel vector in a van der Waals antiferromagnet*, Nat. Commun. **16**, 8911 (2025)](https://doi.org/10.1038/s41467-025-63966-2), published **7 Oct 2025** | Experimental electrical 180° Néel-vector reversal by layer-dependent SOT in **odd-layer A-type (Fe,Co)3GaTe2** with perpendicular magnetic anisotropy | General layer dependence + Néel switching **already experimental** | If a real MnPS3 channel remains, show *source-specific FE paired states*, correct **honeycomb Mn-sublattice** exchange torque and a causal response not equivalent to standard layer SOT. Our old A-type toy in fact resembles the competitor's texture **more than real MnPS3**. |
| **Haque et al.**, [*Terahertz field-induced giant symmetry modulations in a van der Waals antiferromagnet*, Nat. Commun. (2026)](https://doi.org/10.1038/s41467-026-78519-4), published **20 Oct 2026** | **MnPS3** experiments: strong-field THz driven, time-resolved SHG; phonon modes at **1.7 and 4.5 THz**; large coherent lattice response, dynamical reduction of mirror symmetry, possible field-driven magnetic change | Exact materials family and time-domain/symmetry experimental readout | The paper's *phonon-driven nonlinear nonequilibrium* is **not** our hypothetical *linear, equilibrium electronic E→torque susceptibility*. However a claimed electronic angular-momentum mechanism requires excluding its phonon/charge rearrangement contribution and offering an observable that separates them. Do not advertise generic 'THz symmetry in MnPS3' as new. |
| **Ilyas et al.**, [*Coherent terahertz control of metastable magnetization in FePS3*, Nat. Commun. (2026)](https://doi.org/10.1038/s41467-026-78164-x), published **6 Oct 2026** | Coherently drives infrared phonons in FePS3, with nonlinear phonon coupling to Raman mode and induced metastable magnetization | Analogous phosphorous trichalcogenide spin–phonon driven transient magnetism | Confirms need to model/diagnose spin vs lattice angular momentum channels; **not a published observation of our MnPS3 bilayer electronic torque**. |
| **Kang et al.**, [*Octupole-driven spin torque switching of antiferromagnetic tunnel junctions*, *Nature* (2026)](https://doi.org/10.1038/s41586-026-11172-5), published **21 Oct 2026** | Experimental **PtMn3/MgO/PtMn3** AFM junction switching with octupole-mediated spin-transfer torque, modelled by **imbalance of intra- and inter-sublattice spin currents** | Sublattice-resolved angular momentum balance **already produces AFM switching** | Metallic tunnel geometry, noncollinear cluster octupole and dc transport differ from insulating polar **collinear honeycomb** bilayer. Do not claim "inter-/intra-sublattice spin-transfer imbalance drives AFM" as universal novelty. Must show a real material-specific FE switchable **mode-projected** response and proper experimental observable. |

The four articles were already indexed P30–P33 in [the cumulative peer-work audit](PEER_WORK_AUDIT_20261021.md); **this note refines them into operational exclusion criteria**.

## 2. Specific correction to the *source* MnPS3 paper: Eq.13 vs relativistic SOC

The [Sun et al. primary paper](https://doi.org/10.1038/s41467-026-69635-2) writes for its **L2, path 1** pair
\[
[E\parallel M_z M_y]\,P_1^+(L_2)=P_1^-(L_2)\quad\text{(Eq.13)}.
\]
Here the **first component E is an independent nonrelativistic spin-space operation**, while the spatial product \(M_zM_y\) is, as an ordinary orthogonal operation, \(C_{2x}\) in a conventional Cartesian frame with mirrors normal to z and y. Its spatial action is \((x,y,z)\mapsto(x,-y,-z)\), and **physical spin as an axial vector under the corresponding proper rotation** transforms as \((s_x,s_y,s_z)\mapsto(s_x,-s_y,-s_z)\).

**Do not assume that the nonrelativistic spin-space \`E\` is the spinor transformation of a SOC-included electronic state.** With SOC, only an independently established magnetic space-group operation (possibly involving time reversal and a spin rotation) can relate **full spin-operator and torque Kubo tensors** between the \(\pm P\) structures. Its actual site permutation and moment-vector transformation must be checked against the source R1 four-Mn site map. A geometric \(M_zM_y\) check and a nonrelativistic spin-band reversal do not establish that the SOC matrix-valued response simply changes sign at every frequency.

The *source paper* applies Eq.13 to **nonrelativistic spin-band and Boltzmann spin-resolved conductivity classification**. It does not supply, in the accessible main text, a fully validated SOC-covariant electric-field-to-local-Mn-exchange-torque retarded susceptibility. **Absence of that formula from the main text is not evidence that no equivalent calculation exists in its supplement or wider literature.**

An independently executed, **synthetic-only** algebra script is now available: [four-Mn mapped-site operator preflight](../scripts/mnps3_l2_site_operator_gate.py). It checks the labels and exact site bijection and **compares** the nonrelativistic \`spin_E\` action with the axial spinor \`spin_C2x\` action; it **does not decide** which operation is a genuine relativistic symmetry. The script's artificial self-test is **not a measured L2 site map**.

## 3. Decision and falsification criteria

1. **Before any future job:** recover the user's private R1 **four Mn site/moment/structure mapping**, with SHA-256 provenance and correct plus/minus permutation. Passing the synthetic script alone is not an R2 material acceptance.
2. **Build two physically different order-parameter observables:** ordinary layer-odd electronic spin density and **site-weighted Néel-mode electronic generalized force** \(F_{\theta_\nu}=-\partial H_\mathrm{ex}/\partial\theta_\nu\). They are not the same for MnPS3. Only the latter, coupled to a source-identified magnetic mode and dynamical equation, is plausibly an experimental Néel response.
3. **Relativistic symmetry gate:** classify the actual **SOC magnetic group for both source-supported polar stackings**; get allowed field-driven force components without importing Eq.13's nonrelativistic spin-operation as a physical tensor identity.
4. **Spectral/origin gate:** compare purely electronic interband, ionic/phonon, standard SOT/NSOT and transport contributions in a genuine controlled and measurable frequency protocol. The formal operator continuity identity and nearly cancelled toy torques are not unique mechanisms.
5. **If site mapping or operator/gauge exports cannot be validated:** return BLOCKED rather than moving to full Wannier, magnetic damping, or claiming novelty.

**Current platform state:** R1 **CLOSED / NOT NUMERIC_CONVERGED**. This is a **literature and model-operator planning record only**. There is no authorized Slurm/WSL/Japan new job in this note.
