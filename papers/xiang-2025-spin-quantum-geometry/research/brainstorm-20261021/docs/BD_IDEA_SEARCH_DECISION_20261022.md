# B+D scientific idea search: what deserves the next calculation?

**Date: 2026-10-22. Decision memo / literature-based inference, not a confirmed new effect.**  
**Project status:** Source-faithful MnPS3 R1 site/structural checks passed but its two DFT relaxations **did not converge** (zero completed ionic steps); R1 is **CLOSED**. The user's newly approved R2-P0 work is **read-only L2 four-Mn/site/mode source validation**, not VASP restart or compute spending. All prior Gate1–6 outputs are **synthetic**, not MnPS3 predictions.

## 1. The final scientific question, rather than a workflow goal

Can a **physical ferroelectric stacking path**, with the **same controlled magnetic order**, change **which antiferromagnetic *normal mode*** can be electrically excited, in a way that is **not inferable simply from** (i) the spin-band reconstruction, (ii) layer Edelstein/current spin, (iii) ordinary static spin magnetoelectric tensor, (iv) standard collinear-AFM SOT/sublattice torque classification, or (v) phonon-induced THz symmetry changes?

This is a **falsifiable hypothesis only**. A merely FE-even/FE-odd susceptibility, a generic nonzero Néel torque, an SOC-only activation, a spin Hall/NMEE redistribution, or a Ward-identity near-cancellation is **not** a novel paper-level statement.

## 2. Closest prior art — fresh publisher checks added 2026-10-22

| Article (full DOI) | Already establishes | What B+D must NOT claim as new |
| --- | --- | --- |
| [Zhang, Wang et al., *Deterministic Switching of the Néel Vector by Asymmetric Spin Torque*, **PRL 136, 096702** (2026)](https://doi.org/10.1103/fkyr-z5b8) | Sublattice-unequal spin accumulation produces asymmetric field-like/damping-like Néel torque and switching | "Nonidentical Mn-sublattice torque controls Néel order" |
| [Song et al., *Spin-orbit torques in bulk collinear antiferromagnets: Complete classification and induced spin dynamics*, **PRB 113, 224432** (2026)](https://doi.org/10.1103/ktpw-7s2b) | Six MPG SOT classes; Kubo sublattice SOT; Landau–Lifshitz–Gilbert switching | Another generic magnetic-symmetry tensor classification or electric Néel torque |
| [Islam, Ahmad & Xue, *Spin-orbit torque switching of Néel order in band-inverted antiferromagnetic bilayer MnBi2Te4*, **PRB 113, 134440** (2026)](https://doi.org/10.1103/ljqk-968f) | **Insulating in-gap interband SOT** and switching in AFM bilayer | "A gapped 2D AFM can have electric torque without carriers" |
| [He et al., *Chirality and polarization of inertial antiferromagnetic resonances driven by spin-orbit torques*, **PRB 112, 224423** (2025)](https://doi.org/10.1103/1cjr-7cgl) | Mode polarization/handedness changes by adjusting AC SOT polarization | Generic "electrically select a magnon polarization" |
| [Xiang, Jin & Wang, *Third-order intrinsic anomalous Hall effect in altermagnets*, **PRB accepted 19 Oct 2026**](https://doi.org/10.1103/8nkd-ksnd) | A **SOC-enabled, spin-group-forbidden** high-order intrinsic Hall coefficient with Berry-curvature quadrupole; accepted paper rather than an established final issue citation | Simple "spin-group forbids but SOC opens a geometric signal" |
| [Meguro, Ishizuka & Nomura, *Nonlinear spin-motive force driven by mixed-space quantum geometry*, **PRB 114, 194416** (19 Oct 2026)](https://doi.org/10.1103/gdnn-wrbc) | Mixed (k,m) Berry/metric electronic geometry drives magnetization-dynamics-induced AC/DC responses even in an insulator | Mixed magnetization–momentum quantum geometry as a new generic mechanism |
| [Xiang, Jin & Wang, *Spin Transport Revealed by Spin Quantum Geometry*, **PRL 135, 146303** (2025)](https://doi.org/10.1103/14mp-263q) | Spin/Zeeman QGT as framework for SHE/Edelstein; **linear AC displacement spin Hall in insulators**, nonlinear Drude spin responses | Subgap AC spin currents or spin QGT as a new general effect |
| [Sun et al., *A unified symmetry framework for spin–ferroelectric coupling in altermagnetic multiferroics*, **Nat. Commun. 17, 3101** (2026)](https://doi.org/10.1038/s41467-026-69635-2) | Three FE switching classes and **MnPS3** state connections/spin-band/Boltzmann signatures | FE path/symmetry class alone as new |
| [Yang et al., *Symmetry-Programmed Pure Spin Photocurrents in Sliding Ferroelectric Altermagnets*, **ACS Nano** (2026)](https://doi.org/10.1021/acsnano.6c06799) | Sliding FE-paired MnPS3 spin photocurrent tensor | Any FE-programmed spin photocurrent in MnPS3 as untouched |

Other critical references remain in [the cumulative audit](PEER_WORK_AUDIT_20261021.md): Guo vdW layer SOT, Tang insulating AFM NSOT, Zhu stacking magnetism connection, Haque driven MnPS3 THz phonon/symmetry, and Kang octupolar sublattice spin-current imbalance. **Scope of checks:** publisher abstract/main article metadata and previously reviewed materials; absence of an exact mode-projected formula here does NOT prove that no equivalent exists in supplements or other papers.

## 3. Ranked idea decision — high scientific information per unit compute

### A (PRIORITY): FE stacking *selects magnetic eigenmodes*, not just torque sign

**Question:** For the actual **four Mn sites** in a bilayer, do distinct published polarization-switching paths change the **allowed electrodynamic drive of different Néel eigenmodes** (bright/dark, source-selective, or relative spectral oscillator strengths) while leaving a conventional integrated or total-spin signal nondiagnostic? The proposed measurable object is a **mode-projected generalized exchange force**, not the layer-odd toy spin:
\[
F_\nu=-\partial_{\theta_\nu}H_{\rm ex},\qquad
\Xi_{\nu b}(\omega)=\delta\langle F_\nu\rangle(\omega)/E_b(\omega).
\]
A physical mode is defined by a properly normalized eigenvector of **source-validated local-moment exchange + anisotropy dynamics**, not merely by naming a linear combination of electronic spin sites. Its experimental readout would have to specify the mode spectral susceptibility (e.g. field-controlled AFM resonance oscillator strength and selection rules); **a computed \(\Xi_{\nu b}\) without coupling to magnetic modes is only a drive vertex**.

**Lowest-cost falsification BEFORE any DFT:**
1. R2-P0 extracts the *real* four-Mn \(L_2\) site map in both \(P_1^\pm\), and source Eq.13 nonrelativistic vs true **SOC magnetic** connection. Do not assume \(S_{\rm layer-odd}=S_{\rm Neel}\).
2. Use magnetic group stabilizers of each FE daughter, **not** only their pair connection, to enumerate linearly allowed \(E\to F_\nu\) and distinguish genuine zero-to-nonzero mode accessibility from mere partner rotations/sign changes.
3. Include a **second published FE path** solely at source/operator level (not second DFT job), because one partner pair and known connection will often yield only mandatory parity. Are there two pathways with different physical mode-drive matrices not reducible to axis/basis rotation?
4. Test mode branch existence vs known [Song 2026 six SOT classes](https://doi.org/10.1103/ktpw-7s2b), [Zhang 2026 asymmetric SOT](https://doi.org/10.1103/fkyr-z5b8), and [He 2025 resonant SOT chirality](https://doi.org/10.1103/1cjr-7cgl) **at the observable/equation level**.
5. **KILL if** all transformations are simple magnetic-group covariance already implied by earlier classifications, or no mode survives magnetic anisotropy/known even/odd constraints, or no feasible readout exists. Passing a toy site-bijection script is **not** sufficient.

**Physical upside:** FE switching controls **which dynamical AFM degree of freedom couples to AC E**, not merely reversing the sign of an existing spin current. **This remains a plausible but extremely competitive research angle**, not an established gap in literature.

### B (BACKUP): layer/actual-site-resolved *insulating spin displacement* and proper angular-momentum continuity

**Question:** Does an experimentally addressable finite-frequency **interlayer spin/angular-momentum transfer**, with full electronic exchange + SOC torque and physical site boundaries, reveal a FE-state spectral difference that cannot be reduced to Xiang's **linear displacement SHE**, conventional Néel SOT, or changing a model-specific spin-current partition?

**Low-cost gate:** derive a gauge-invariant proper observable with an explicit boundary/torque-dipole contact term, compare spectral sum rules and \(t_\perp\to0\), Fermi-gap, \(\eta\to0\), and Wannier rephasing. If only the Ward identity or a compensated torque remains, **KILL**. Do not start from \`wannier90_hr.dat\` alone.

**Risk:** physics definition and experimental readout are difficult, and closest competitors already contain AC insulator spin transport and Néel torques. Treat as lower-priority exploratory theory, **not** a reason to start new DFT.

### C (SUPPORTING MECHANISM, not stand-alone paper): spin-space-to-relativistic symmetry reduction

Actual source **L2, path1** [E || MzMy] is a nonrelativistic spin-space relation; physical \(M_zM_y=C_{2x}\) rotates an axial spin as \((s_x,s_y,s_z)\mapsto(s_x,-s_y,-s_z)\). With SOC there may be a distinct antiunitary or no exact daughter connection, depending on **source four-Mn pattern and moment axis**. Establish this carefully, because a spurious switching parity destroys all later physics.

**Low-cost gate:** site-dependent spinor/magnetic double-group analysis with \(L_2\) and \(P_1^\pm\); compute allowed tensor components **for the true MPG**, compare leading SOC order only if formally meaningful and the symmetry permits it.

**KILL as novelty if** it is just standard SOC reduction; [OSSG Nature 2026](https://doi.org/10.1038/s41586-026-10401-1) and [Xiang et al. accepted PRB 2026](https://doi.org/10.1103/8nkd-ksnd) already occupy the general spin-group-vs-SOC concept.

### D (DO NOT INVEST NOW): broad new MnPS3 DFT/NEB/Wannier or more arbitrary toy Hamiltonians

R1 reached Slurm timeout without a single converged electronic step. A major multi-path VASP calculation now would **spend resources without a sharp observable**. The prior 32-band synthetic A-type AFM differs from MnPS3 site Néel order; its compensated torque has no verified MnPS3 dynamical readout. **No new jobs** unless Idea A passes a symmetry/mode-selection and literature criterion; if it does, start a NEW user-approved bounded SCF stabilization campaign, not rerun historical R1.

## 4. Single next scientific deliverable (read-only)

Produce an **MnPS3 L2 FE-path magnetic normal-mode selection matrix**:

| FE state/path | four Mn site sign/order from real input | SOC MPG + stabilizer | local moment normal-mode irrep / trial mode (not yet frequency) | \(E_x,E_y\) coupling allowed? | actual equality/zero relation across FE switch | closest SOT/resonance prior-art formula | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| \(P_1^+(L_2)\) | **requires private source** | ? | ? | ? | ? | Song/Zhang/He | BLOCKED until verified |
| \(P_1^-(L_2)\) | **requires private source** | ? | ? | ? | ? | Song/Zhang/He | BLOCKED until verified |
| Path2 comparator | published, no new DFT input | ? | ? | ? | ? | Song/Zhang/He | SOURCE-PREFLIGHT ONLY |

This is a **scientific decision table, not an invented material property**. If the group algebra cannot distinguish a real physical mode, report **NO-GO for A** and pivot to a less saturated topic rather than reflexively designing an expensive calculation.

**Research narrative test:** "Stacking FE acts as a *mode-addressability switch* for compensated antiferromagnetism, with source-specific null controls and a measurable resonance." That is a *candidate story*. It only becomes good physics if the local-moment eigenmode, causality, competing papers, and actual material signal can be independently established.

**Infrastructure status:** [tasks/CURRENT.md](../tasks/CURRENT.md) authorizes ONLY the **read-only R2-P0 site/magnetic-operator audit**. Writing this memo does not start DeepSeek or HPC. R1 remains CLOSED; no R2 compute permitted.
