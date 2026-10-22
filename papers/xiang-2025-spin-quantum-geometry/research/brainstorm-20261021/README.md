# Independent A/B/C/D brainstorming (2026-10-21)

## R1 real-material terminal status (2026-10-22)

**R1 CLOSED BY USER; no active B+D material jobs or authorized restart.** The published MnPS3 L2 path-1 geometry/magnetic mapping preflight passed, but **both real paired Slurm relaxations reached two-hour TIMEOUT in their first electronic SCF**; no ionic relaxation, converged energies, FE property, layer spin or torque have been obtained. `RUN_FINISHED=true` means only that the scheduled jobs have ended, while `NUMERIC_CONVERGED=false` and `FE_VERIFIED=false`. The user explicitly accepted this failed numerical result as the R1 endpoint. **Do not auto-resubmit from a cached historical command.**

[Sanitized scientific closeout](docs/BD_MNPS3_R1_CLOSEOUT_20261022.md) · [Live task pointer: CLOSED / no active job](tasks/CURRENT.md) · [Historical R1 instruction, do not execute](tasks/BD-MNPS3-L2-P1-EXECUTE-R1.md). Private signed run evidence remains in the user's private research archive. Any R2 restart requires its own approved task and convergence-control plan.

**RESEARCH EXPLORATION ONLY — NOT Xiang–Jin–Wang PRL 135, 146303 (2025) reproduction.**
This folder archives user/assistant explorations motivated by that paper. The original `paper_models.py`, `methods/`, `tests/`, frozen outputs and scientific acceptance state remain unchanged.

## Distinct scientific statuses

- **A:** parked / falsifiable. Nonlinear hidden layer-spin response must pass continuity, finite-q and boundary-observable gates.
- **B:** toy-model checks pass for the phase-only null test and physical-parameter modulation; no novelty claim. The zero BZ integral in the null test holds only under its smooth periodic, k-independent conjugation hypotheses.
- **C2:** four-band AC interlayer spin-transfer response and Ward identity verified; *NO-GO for the original novelty claim* because conventional mixed Berry curvature and prior insulating Néel spin-orbit torque account for the effect.
- **D:** competing switching-operator constraints are Neumann/stabilizer algebra, not a new theorem. Physical magnetic connection and phase-resolved observables remain to be tested.

## Research workflow: three separate records (2026-10-21)

1. **Idea:** A/B/C/D hypotheses are exploratory; read [ABCD round 1](docs/ABCD_round1.md). No new mechanism is established by proposing a tensor or connection operator.
2. **Testing:** [Registry-built PT magnetic bilayer Gate 1](docs/STACKING_REGISTRY_GATE1_20261021.md), [script](scripts/stacking_registry_gate1.py). This synthetic first step checks physical-bond interlayer hopping, orbital-embedding seam, PT/C2z partner relation, spin-torque identity, metallic/insulating controls and the necessary `t_perp=0` falsification. **NO-GO for the first interlayer novelty claim**; a switchable Drude proxy mainly survives as a disconnected-layer effect in the optional assumed chiral channel.
3. **Peers:** [Verified peer-work / competing-paper audit](docs/PEER_WORK_AUDIT_20261021.md) records direct overlaps, DOI sources, uncertainty levels, and the remaining gates. Do not turn abstract-level overlap into claimed equation-level proof.

**Scientific boundary:** nothing here changes `paper_models.py`, `methods/`, paper `tests/`, paper `results/`, paper `docs/ACCEPTANCE.md` or original-prl reproduction scientific status. This directory is **not published-paper data**. `C2z` registry partners are not established sliding-ferroelectric states.

## Documents

- [A exploration and counterexamples](docs/A_exploration.md)
- [B/C deep assessment](docs/BC_research_assessment.md)
- [ABC → C2 scientific handoff](docs/ABC_to_C2_handoff.md)
- [C2 gate-0/2 NO-GO audit](docs/C2_no_go.md)
- [A/B/C/D round 1](docs/ABCD_round1.md)

## Independent scripts (toy models, not original paper reproductions)

- [C2 minimal test](scripts/c2_minimal.py): finite-frequency Kubo/Liouville, interlayer/local torque Ward identity, disconnected limit.
- [B rephasing null test](scripts/b_geometry_nulltest.py) and [three-parameter Bianchi check](scripts/bc_three_parameter_bianchi.py) are independently runnable; they are not paper reproduction tests.
- [PT doublet projector and Drude reality check](scripts/pt_hidden_metric_and_drude.py): new conditional PT-staggered interband kernel and Fermi-surface model; no novelty claim.

## Round 2: PT-staggered spin-geometry audit

- [Round 2 report](docs/PT_staggered_spin_geometry_round2.md): U(2)-invariant doublet projector calculation; distinction between a surviving geometric kernel and a vanishing insulating dc Drude response. **Exploratory; novelty not established.**

## Research integrity and provenance

Numerical toy-model tests do **not** establish a new quantum-geometric effect or any material-specific realization. Original papers and supplemental PDFs are not mirrored. Publication-priority screening is provisional; check linked literature before claiming originality.

This directory deliberately has **no imports into paper reproduction code and no new CI acceptance claim**.


**2026-10-21 Gate 1b correction:** [Complete four-band PT operator-basis enumeration](scripts/pt_operator_basis_gate.py) leaves only six PT-even Hermitian `tau_a sigma_b` products. Direct spinful interlayer `tau_{x,y} sigma_i` terms are forbidden **in this specific four-band PT representation**. The next scientifically meaningful construction requires extra orbital freedom or justified PT breaking; see the amended [Gate 1 report](docs/STACKING_REGISTRY_GATE1_20261021.md).

## B+D Gate 2 (2026-10-21): eight-band physical-registry experiment

- [Scientific report](docs/BD_EIGHTBAND_GATE2_20261021.md) and [standalone NumPy code](scripts/bd_eightband_gate2.py).
- Two opposite inversion-parity orbitals per layer permit PT-even spin-dependent interlayer hopping; every hopping term is generated from displaced real-space bonds. This is **symmetry-compatible toy hopping**, not first-principles Wannier/Slater–Koster extraction.
- With both ordinary and spin-dependent hopping present, the gapped model has an interband AC layer-staggered **spin-density** response. It changes sign across C2z-related registries and vanishes when either hopping family is absent, unlike the Gate 1 disconnected-layer artifact.
- Independent rank-two PT-doublet projectors, random U(2) mixing, orbital-embedding k sewing, direct/indirect gaps, k convergence and torque commutators tested.
- **NOT an FE pair, NOT a new quantum metric, NOT peer-priority clearance, NOT a completed original-paper reproduction.** Next gates: microscopic orbital SOC hopping, magnetic connection to actual polar states, torque/edge observable, equation-level peers.


## B+D Gate 3 (2026-10-21): atomic-SOC microscopic path and negative novelty tests

- [Scientific source/peer/physics report](docs/BD_GATE3_ATOMIC_SOC_MAGNETIC_NO_GO_20261021.md) and [independent standalone 16-band NumPy script](scripts/bd_gate3_atomic_soc.py).
- **PASS (model):** explicit spin-INDEPENDENT \(s,p_x,p_y,p_z\) Slater–Koster interlayer bonds plus onsite atomic \(L\cdot S\) produce PT-even spin-dependent \(s_T\to p_{z,B}\) hopping after \(p_x,p_y\) Schur downfolding. Controls \(V_{sp}=0\) and SOC=0 remove it. This is **well-known perturbation theory**, not novelty.
- **NO-GO (original magnetic claim):** subgap AC layer-staggered \(S_-^z\) response persists at \(m=0\) and is even under Néel reversal.
- **OVERLAPS (geometry):** static \(E_x\to S_{\rm tot}^z\) signal is magnetic-odd, but independently agrees with **standard mixed \((k_x,b_z)\) Berry curvature** (difference ~\(1.3\times10^{-12}\) on \(11^2\)).
- **SYMMETRY TEMPLATE, NOT FE:** exact PT forbids *net* \(P_z\). An \(m_zT\) *connection* maps two PT-broken daughter models with manually imposed layer potentials \((+\delta,+u)\leftrightarrow(-\delta,-u)\); total and staggered spin response have opposite switching parities. No self-consistent polarization, structural double well or FE barrier demonstrated.
- Next decisive gate: genuinely polar sublattice geometry + hopping/charge-energy-derived polar asymmetry, **not an inserted \(u(\delta)\)**; then full proper spin/torque/finite-\(\omega\) response and near-neighbor equation-level comparison.

The original PRL-paper reproduction, frozen data, CI and acceptance statuses remain strictly untouched.


## B+D Gate 4 (2026-10-21): genuine A/B registry and symmetry-filtered spin response

- [Detailed science / peer review](docs/BD_GATE4_POLAR_HONEYCOMB_MAGNETIC_RESPONSE_20261021.md); [structural/electronic polar honeycomb toy](scripts/bd_gate4_polar_honeycomb.py); [32-band orbital/SOC/AFM spin Kubo](scripts/bd_gate4b_atomic_spin.py).
- **Structural test:** inequivalent A/B honeycomb sublattices, real displaced bonds and an *uncalibrated* band+screened-pair energy produce two opposite-electronic-dipole AB/BA **local minima** and a candidate saddle. Charge transfer vanishes for disconnected layers or identical A/B species. Ionic-model parameters matter.
- **Physical group constraint:** stable polar AB/BA retain \(C_{3z}\), which **strictly forbids linear** \(\chi_{S_z,E_x}\) and \(\chi_{S_z,E_y}\), despite generic-registry prior nonzero toy signals.
- **Magnetic readout model:** 32-band all-spinless spatial Slater–Koster hopping plus onsite atomic SOC and fixed A-type AFM (hypothetical material). An **\(m_zT\) connection** gives FE-even total \(\chi_{S^x,E_x}\) and FE-odd layer-staggered \(\chi_{S_-^x,E_x}\). The response vanishes for magnetic \(m=0\), atomic SOC=0 or interlayer hopping=0. **Not novel standard symmetry or magnetoelectric physics.**
- **Important negative result:** in the 32-band + empirical ionic toy energy, AB/BA are only metastable; AA is lower. No physical switching barrier/real FE claim.
- **Peer-work correction:** 2026 *Nature Communications* already has a three-type FE-spin symmetry framework with **six DFT FE minima in bilayer MnPS3**. 2022 npj already has FE-controlled light-induced spin/orbital magnetization. See [expanded peer audit](docs/PEER_WORK_AUDIT_20261021.md).

All files remain under independent \`research/\`. Xiang et al. PRL reproduction frozen code, tests, CI and scientific acceptance are unchanged.


## B+D Gate 5 and compute handoff (2026-10-21)

- [Subgap interlayer spin transfer vs local-torque Ward audit](docs/BD_GATE5_TORQUE_WARD_20261021.md) and [independent NumPy test](scripts/bd_gate5_torque_ward.py). Synthetic 32-band AFM polar registry: at \(\omega=.06\), the imaginary interlayer-transfer susceptibility is about \(+1.02174\times10^{-4}\), while local torque is \(-1.02756\times10^{-4}\), net only \(-5.81921\times10^{-7}\) (23×23 mesh). The exact response continuity residual is ~\(10^{-17}\); opposite polar partner changes signs; interlayer/SOC/AFM-off null tests pass. This is a **standard Ward/torque-cancellation model result, not a novel material effect**.
- [GPT / GitHub / local DeepSeek harness / WSL-HPC scientific responsibilities and handoff](docs/BD_RESEARCH_TO_COMPUTE_HANDOFF_20261021.md). Uses the existing group theory, structure core/benchmark, BSF generator, private Wannier/NLO workflows, private evidence archive and local bridge. Distinguishes \`IDEA\` / \`MODEL_PASS\` / \`JOB_FINISHED\` / \`NUMERIC_CONVERGED\` / \`PROPERTY_ACCEPTED\` / \`NOVELTY_REVIEWED\` with SHA evidence. Real WSL/Japan jobs and spin-response acceptance are **not** claimed to have been executed in this conversation.
- GitHub is the **auditable research record**; raw licensed/private DFT data and scheduler credentials stay on private compute storage. The manuscript reproduction's frozen scientific status is unchanged. Review/test the actual independent spin/layer/position matrix elements before reusing a Wannier \`hr.dat\` for spin-geometry physics.


## Real-material R0: MnPS3 bilayer (2026-10-21)

- [Real-material source/preflight and precise DeepSeek handoff](docs/BD_MNPS3_REAL_MATERIAL_PREFLIGHT_20261021.md). **NOT SUBMITTED**; no project DFT, WSL/Japan data or confirmed source POSCAR pair yet.
- Primary source: Sun et al., *Nat. Commun.* **17**, 3101 (2026), [DOI](https://doi.org/10.1038/s41467-026-69635-2). **Chosen first pair:** published `P_1^+(L_2) ↔ P_1^-(L_2)` (path 1, **author Eq.13** `[E || M_z M_y]`), **not** the toy's arbitrary `m_z T`. Authors state L2 ground-state magnetic ordering and 11.8 meV/f.u. path-1 barrier.
- Specific novel-physics **hypothesis**, not assertion: material-specific full causal *layer spin density ↔ electronic interlayer spin transfer ↔ local exchange/SOC torque* frequency response, including exact Ward/connection relations and spin-vertex same-gauge completeness. The published spin-band, nonlinear Hall and pure spin-photocurrent work is a severe novelty boundary.
- Preflight task `BD-MNPS3-L2-P1-STRUCTURE-R0` is **structure acquisition / exact moment and symmetry validation only**. WSL/Japan jobs need the user's later explicit approval; never infer that the DeepSeek bridge is connected from its README.

- [Copy-ready local DeepSeek preflight task (R0; **NOT SUBMITTED**)](tasks/BD-MNPS3-L2-P1-STRUCTURE-R0.md) freezes source Eq.13, candidate states, exact layer/moment checks, SHA evidence and no-VASP/no-HPC permissions.


## DeepSeek single-link activation (2026-10-21)

**Stable task URL:** [tasks/CURRENT.md](tasks/CURRENT.md) → [versioned R1 two-job execution spec](tasks/BD-MNPS3-L2-P1-EXECUTE-R1.md).

This is the single GitHub link to pass to the already-running local DeepSeek harness. The R1 mission obtains original MnPS3 L2 path-1 structures, checks exact magnetism and physical source, performs signed dry-runs and, **if all structural + platform + resource/approval gates pass**, submits only two bounded matched WSL/Japan baseline runs with explicit job IDs, provenance and native convergence audit. No automatic Wannier/NLO run, no invented HPC allocations and no confidential output in this public tree.

**Status on GitHub: READY FOR LOCAL AGENT INTAKE / NOT RECEIVED / NOT SUBMITTED.** The agent or live queue must report actual execution status. Updating files is not submission.


## Gate 6: exchange vs SOC torque — 2026-10-22

- [New independently tested operator decomposition](docs/BD_GATE6_EXCHANGE_SOC_MAGNETIC_TEXTURE_AUDIT_20261022.md), [standalone NumPy verification](scripts/bd_gate6_exchange_soc_torque.py) (imports original Gate4 toy files). In the same toy AB registry at ω=.06 (23×23), **interlayer electronic spin-transfer ≈ +i1.021742e-4**, **exchange electronic torque ≈ −i1.008285e-4**, **atomic SOC torque ≈ −i1.927658e-6** and **net ≈ −i5.819214e-7**. The operator Ward identity holds to ~7e-18; no claim that an actual Mn magnetic mode is driven.
- **Critical physical correction:** Gate4/5 synthetic AFM is *A-type layer-staggered* \`m tau_z sigma_z\`; **real MnPS3 has two opposite-spin Mn honeycomb sublattices within each layer**, so layer-odd spin is not the physical L2 Néel operator. Define source-mapped four-Mn sublattice projectors and \`F_theta=-dH_ex/dtheta\` before promoting any predicted layer torque to a magnetic-order response.
- **Competition tightened:** 2025 *Nature Communications* layer-dependent vdW SOT; 2026 *Nature Communications* **20 Oct** THz-driven MnPS3 phonon/symmetry response; 2026 *Nature* **21 Oct** sublattice/octupole spin-transfer-torque switching in an AFM tunnel junction. Interlayer/local torque compensation, THz dynamics and spin-current imbalance are **not unoccupied general mechanisms**.
- [Zero-compute read-only MnPS3 L2 site/mode contract and R1 SCF diagnostic plan](tasks/BD-MNPS3-R2P0-SITE-OPERATORS-DRAFT.md), **DRAFT ONLY**, not a job command. **R1 remains CLOSED**, and no R2 computing or restarts are approved/submitted.


## R2-P0 approved read-only review (2026-10-22)

The user approved **proceeding with read-only scientific analysis** of existing R1 evidence, **NOT new VASP/Slurm/WSL/Japan/para computation**. The stable [CURRENT task pointer](tasks/CURRENT.md) now links the [actionable R2-P0 instructions](tasks/BD-MNPS3-L2-R2P0-READONLY.md). The **local DeepSeek harness has not yet provided an execution receipt**.

- [Publisher-verified four-paper 2025–26 competitor/novelty review](docs/BD_PEER_UPDATE_20261022_SOT_THz_OCTUPOLE.md): layer SOT in A-type FCGT; THz phonon-driven MnPS3 symmetry; coherent FePS3 magnetization; octupolar PtMn3 spin-torque tunnel junction.
- [Four Mn site-projector algebra, source-hash gating and spin-space vs axial-spin checks](scripts/mnps3_l2_site_operator_gate.py): synthetic self-test performed; **no source-correct four-Mn map has yet been returned**, and the simple matrix representation is **not** accepted SOC magnetic group or PAW/Wannier spin vertex.
- [Gate 6 physics analysis](docs/BD_GATE6_EXCHANGE_SOC_MAGNETIC_TEXTURE_AUDIT_20261022.md) remains model-only; unlike that *layer-A-type* model, actual MnPS3 has intralayer honeycomb Néel sublattices. The R2-P0 task is to recover R1 private plus/minus mapping, define correct \`F_theta=-dH_ex/dtheta\`, screen symmetry components and propose controlled SCF stabilization **without executing it**.
- **R1_CLOSED remains immutable:** two real jobs ended at 2-hour scheduler TIMEOUT in the first SCF, neither electronically converged, no accepted magnetic moment/FE or response output.


## Scientific reprioritization: 2026-10-22 idea-search decision

[**Ranked idea and closest-new-publication kill-test memo**](docs/BD_IDEA_SEARCH_DECISION_20261022.md). 2026 full SOT classifications, asymmetric Néel torques, insulating AFM bilayer in-gap torques, resonance handedness control, new Xiang SOC-only intrinsic Hall, and mixed k/m geometry **remove generic novelty** from earlier B+D proposals. The only currently highest-ROI **candidate** is a *source-verified FE switching path × actual MnPS3 four-Mn normal-mode* **electric driving selection matrix**, **if it is not reducible** to those established works. The alternative proper layer spin displacement/torque observable remains backup with difficult experimental interpretation. The next useful artifact is a symmetry/mode-based \`ALLOW/FORBID\` decision table with a source-correct L2 site map, not another unfitted model or DFT. No new job submission authorized; R1 remains CLOSED, R2-P0 remains read-only.


## 2026-10-22: Full layer-group theoretical pivot, peer-first

**[General-theory hypothesis / exact unitary-domain no-go and stop gates](docs/BD_GROUP_THEORY_80_528_PIVOT_20261022.md).** The group-theory prototype has been uploaded to the existing \`group-theory-operations-toolkit\` in [protected-main Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32), NOT merged or advertised as a complete 80/528 catalog. A standalone NumPy test verifies that the allowed real static electric-field-to-magnetic-mode coupling for **exact full magnetic unitary symmetry partners** is merely conjugated by the mapped representations, so **rank, singular values and allowed intertwiner-space dimension cannot be created by a FE-domain coordinate rotation**. A lab-fixed off/on may be entirely orientation-only.

New direct competitors include the [2026 complete spin-layer-group corepresentation PREPRINT](https://arxiv.org/html/2605.25484v2), [FE switched magnon transport, Nano Lett 2026](https://doi.org/10.1021/acs.nanolett.6c02164), and [PT electric chiral magnon MLG PREPRINT](https://arxiv.org/abs/2601.16052). The only honest next target is to distinguish exact-magnetic-domain covariance from (i) **SOC-obstructed magnetic connections**, (ii) **nonconjugate magnetic/site stabilizers**, and (iii) **finite-q nonsymmorphic site/mode selection with a real momentum source**. No new material result, all-MLG census or independent novelty has been established; no new VASP jobs authorized.


## 2026-10-22 overall review: new antiunitary and peer-competition gate

**[Cross-project scientific audit](docs/BD_CROSS_PROJECT_DEEP_AUDIT_20261022.md)** evaluates the four independent paper-reproduction \`PARTIAL\` states, Gate1–6 toy physics, MnPS3 R1 nonconvergence, 2025–2026 peer competition, and the next publishability/no-go decisions. New direct peers: [2026 PRL magnon geometry by light](https://doi.org/10.1103/12cl-b9jj), [2026 PRB twisted-bilayer MoS2 SHG symmetry](https://doi.org/10.1103/ldws-7yvc), [2026 npj Quantum Materials experimental magnetic vs spin-group phonon selection](https://doi.org/10.1038/s41535-026-00857-9).

**Scientific update to protected-main [group toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32):** a full 2-/4-band finite-frequency Lehmann test shows that an **antiunitary** full-domain connection gives **Onsager-swapped \`chi_AB(partner)=eta_A eta_B chi_BA(source)\`**, NOT naive same-forward-sign parity. The correction is classic Onsager reciprocity, not new theory, but is a decisive rejection criterion for future FE response switching claims. New script [antiunitary Kubo gate on PR branch](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/antiunitary_kubo_onsager_gate.py).

**Priority:** source-real physical site/mode/projector + true magnetic connection + full causal measurement/source pair and **only then** group-wide response screening; moire local-stacking response textures are a cautiously ranked alternative that must beat existing polar merons and twist SHG. No full LG/MLG census, no verified new effect, no fresh HPC authorization; MnPS3 R1 remains CLOSED, R2-P0 read-only.


## 2026-10-22 G2 | p3 Wyckoff-site and physical-local-response screening

[**Independent G2 science report**](docs/BD_G2_P3_SITE_ORBIT_FROBENIUS_20261022.md) and protected-main [**group-theory Draft PR #32**](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) now contain tested source-based \`p3\` (#65 layer group) **site orbit** representation, a **Frobenius reciprocity independent oracle**, and strong negative controls.

- **Same global C3, different site orbit:** per-layer generic 3d \`E_xy → local axial Sz\` Hom dimension **2**, three special 1a/1b/1c **0**; transverse spin components **4 vs 6**. With two layers counts double. All net and layer-summed axial \`Sz\` signals are forbidden. This is **known induced-representation/hidden local multipole physics**, not new effect; actual Kubo coefficients may vanish.
- **More important rejection constraints:** Hom coefficient count is not mode rank (\`rank K≤2\`); for any static real \`C3\`-equivariant \`K\`, \`KᵀK=λI\`, so total normalized local-response norm is **linear-polarization isotropic**. Complex AC \`K†K\` may permit circular-selective response algebraically, not proven optical/magnon physics.
- **No fake FE switching:** a full parent symmetry must preserve site orbit stabilizer order, so generic \`3d\` vs special \`1a⊕1b⊕1c\` are **not** a proven symmetry-equivalent ±FE pair. Any actual switching between them would need a substantially different structure/energy story. 2018 magnetic multipoles, 2020 hidden effects, 2025 electric trimer moments and 2026 hidden spin polarization predate generic site-hidden-spin claims.
- **Next focus G3:** search source-backed genuine polar magnetic stacking partner operations with ordered site/mode representations; compare exact-unitary, antiunitary Onsager and non-conjugate group cases before exhaustive 80/528 sweeps. No new DFT, Wannier or HPC tasks; MnPS3 R1 remains CLOSED.


## Science Idea Kill Gates and genuinely source-backed G3 (2026-10-22)

- **[8-stage scientific falsification / 20-claim ledger](docs/BD_IDEA_ELIMINATION_GATES_20261022.md)**, [machine-readable 20-entry idea registry](data/BD_IDEA_REGISTRY_20261022.json). Claims are recorded as NO_GO_CLAIM, BLOCKED_EVIDENCE, BLOCKED_NOVELTY or OPEN_HYPOTHESIS; valid METHOD/MODEL tests are preserved. A NO_GO rejects a **novelty claim**, not the field or the software.
- [Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) adds a **research-only Python standard-library idea_elimination_gate.py** and two sample JSON fixtures. It is **unmerged**. Self-test mode is supplied, but no execution of the exact GitHub script is certified in this conversation; programmatic statuses cannot validate user-authored scientific evidence.
- [Next G3 source/field/magnetic-mode hypothesis and source-specific kill tests](docs/BD_G3_NEXT_HYPOTHESIS_20261022.md). Fundamental correction: source-backed sliding FE switching typically involves domain walls; **end-point Seitz operation is not a driven switching trajectory**. The potential magnetic-site/local generalized-force response of a **moving FE wall** must beat 2024 nonrelativistic local Edelstein torque, 2025 PRX local magnetic multipoles/surfaces, 2026 PRL polar-bilayer E-torque and known wall kinetics before any new physics claim.
- **Strong preemption:** 2025 npj *shift photocurrent vortices from polar merons* and 2020 PRB sliding moiré topology/charge pumping already use exactly the broad Berry/stacking narrative. They are registered as NO_GO_CLAIM for **generic versions**; genuinely different observables require a new ID and formula-level audit.
- User-approved MnPS3 R1 CLOSED (two SCF timeouts), separate R2-P0 remains READ-ONLY. **No authorized DFT/Slurm/WSL/Japan restart** from these research files.


## 2026-10-22 G3 science ROI reset: avoid domain-wall DFT

[**Cost/physics decision: no magnetic FE domain-wall DFT yet; symmetry-constrained registry-gradient search**](docs/BD_G3_DOMAIN_DFT_COST_AND_REGISTRY_GRADIENT_20261022.md). Primitive-cell MnPS3 R1 failed 2h SCF; a magnetic wall supercell and full localized torque response are higher-risk. 2026 PRB already computes ferroelectric wall Edelstein, 2024 bilayer CrI3 magnetic FE spin/orbital Edelstein, 2026 PRB slow-phonon Edelstein, 2026 *Nature Physics* WTe2 coherent interlayer-sliding dynamics, and 2026 *Nature Communications* MoSe2 phonon-beating transient magnetism. None supports generic originality in combining those keywords.

**Actionable method-only G3.1** in protected-main [group toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32): [registry_gradient_gate.py](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/registry_gradient_gate.py) checks \`K0 in Hom(E,O)\`, \`dK/dδ in Hom(U⊗E,O)\` on an explicit mirror sign representation. This is known tensor-product symmetry and **only a method gate**. Next scientific test: two true BSF source-backed registry/layer-exchange cases with physical U representation, then exact peer-observable test. Full DW DFT remains deferred; no new compute authorized.

[**Idea register**](data/BD_IDEA_REGISTRY_20261022.json) now 24 claims, preserving new FE-wall spin, magnetic sliding spin/orbital and phonon Edelstein overlap as NO-GO to broad originality, while response-gradient as BLOCKED_NOVELTY.


## 2026-10-22: Research-core reset — physics before narrative architecture

[**What we kept circling, direct-competitor vs Introduction role, and one minimal growth-ready research seed**](docs/BD_CORE_PROBLEM_REFERENCE_ROLES_AND_SMALL_SEED_20261022.md). Key user correction: an idea can grow from one surprising, independently measured response; do **not** require a complete 80/528-LG paper structure before the first result. Classify papers *relative to one specific claim*: direct competitor (same observable/source/protocol and overlap), potential direct requiring equation check, or Introduction context. Do **not** classify all shared-keyword papers as claim-killing prior art.

**A first symmetry-only seed, not novel physics:** nonmagnetic current-driven layer-odd \`Sz_T-Sz_B\` vs relative in-plane slide \`u\`. For an IDEALIZED \`C3z + vertical mx + layer-swapping horizontal mz\` group, the zero-shift spin-current response is forbidden and a unique first-order \`u_x*j_y-u_y*j_x\` coefficient is allowed. Explicit SVD group basis gives one invariant; physically correct layer swap matters. No nonzero electronic coefficient, material FE stability or transport protocol is established. [Layer-exchange representation code](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/physical_layer_exchange_gate.py) and [unique bilinear invariant test](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/stacking_first_seed_gate.py) live in protected-main Draft PR #32.

**Next smallest scientific action:** use the original 2026 LEE PRL *minimal model* to see whether the same response derivative is already present; then one symmetry-consistent SOC/transport calculation of \`chi_{layer-odd Sz,j_y}(u_x)\` with odd-u, mirror, fixed-current vs fixed-field and SOC-off controls. **Not yet executed, no new DFT**. The [machine-readable registry](data/BD_IDEA_REGISTRY_20261022.json) now explicitly scopes NO_GO to individual claims and adds only one new seed (BD25), never grants global topic rejection.


## 2026-10-22 BD25 source-model audit: the proposed spin-z response is not in the LEE minimal Hamiltonian

[**Exact source-model no-go and baseline Boltzmann test**](docs/BD25_LEE_MODEL_AUDIT_20261022.md). Zhou et al.'s published minimal Eq. (3) conserves an in-plane spin-helicity operator, so every band has \(\langle\tau_z\sigma_z\rangle=0\); registry-dependent scalar \(t_\perp(u)\) or \(\Delta(u)\) cannot produce the proposed layer-odd \(S_z\) response. This is **not** a real-material no-go: the correct group-theory invariant remains allowed, but a physical spin-z SOC/valley term and microscopic \(\partial_u H\) are missing. A dimensionless Eq. (2)-style baseline gives \(\chi_{S_-^y,j_x}\) even in \(\Delta\), global \(\chi_{S_+^y,j_x}\) odd, and target \(S_z=0\), with parity residuals \(2.3\times10^{-16}\). Code is in the group-toolkit draft branch. No DFT/HPC.


**Source-backed next model candidate:** the [BD25 model audit](docs/BD25_LEE_MODEL_AUDIT_20261022.md) now identifies the 2020 MoS\(_2\) Slater–Koster tight-binding model as a possible source for registry-dependent interlayer hopping, with a separately sourced SOC extension. It has **not** been implemented; next is only a two-layer \(u=0,\pm\delta\) Hamiltonian covariance test. The generic claim “sliding tunes Edelstein response” is pre-empted by PRB 114, 125129 (2026) and magnetic CrI\(_3\) sliding work, but the narrower nonmagnetic layer-odd \(S_z\) derivative is not yet shown to be identical. Registry updated to 26 claim-specific entries.


## 2026-10-22 research correction: examples do not pre-empt a general theory

[**Abstraction-level prior-art policy and a candidate stacking-response-field theory**](docs/BD_THEORY_LEVELS_AND_STACKING_RESPONSE_FIELD_20261022.md). Prior work is assessed along two axes: scope (S0 realization, S1 mechanism, S2 family classification, S3 general principle) and contribution type. A concrete sliding-ferroelectricity or specific ferroic-order paper can be foundational evidence for a later general theory without pre-empting it; only an exact proposition at the same scope or a broader theorem that implies it pre-empts the general claim.

**New S2–S3 hypothesis BD27:** treat a specified layer-resolved response as a field over the stacking-registry torus. Under a layer-swapping horizontal-mirror Seitz operation, a layer-odd \(S_z\) response to in-plane current obeys \(\chi(-\delta)=-\chi(\delta)\). Candidate corollary: response zeros at fixed registries \(2\delta\in\Lambda\) and no even local Taylor orders. This follows from standard covariance and is **not yet a novel theorem**; exact LEE SM / configuration-space prior-art comparison and verification of each actual Seitz operation are required. The [27-entry registry](data/BD_IDEA_REGISTRY_20261022.json) now records scope level and contribution type.


**Source-checked precedents:** [Ji et al., PRL 130, 146801 (2023), General Theory for Bilayer Stacking Ferroelectricity](https://doi.org/10.1103/PhysRevLett.130.146801) explicitly generalizes earlier sliding-ferroelectric examples to all 80 layer groups. [Pang & He, PRL 135, 116402 (2025), Generalized Neumann's Principle](https://doi.org/10.1103/trhd-kxm1) generalizes FQFE and conventional ferroelectricity over 32 point groups. These support the scope distinction and also raise BD27's novelty bar: our candidate must be a response-specific, nontrivial global constraint, not a renamed general Neumann principle.


## 2026-10-22 BD27/BD28: From examples to a general theory of registry-space response zeros and indices

[**Theory/model report: symmetry-forced zeros, response operator parity and torus winding constraints**](docs/BD27_REGISTRY_TORUS_ZERO_INDEX_GATE_20261022.md). This is the first **cross-family** candidate driven by the user's important distinction: a single known material effect does not pre-empt an S2–S3 general theory. The direct prior frameworks are [2023 BSF](https://doi.org/10.1103/PhysRevLett.130.146801) for six hexagonal G/A/B/C/M/N registry points and [2026 LEE](https://doi.org/10.1103/qhvs-pfkq) for layer-resolved current-induced spin/80-LG selection. The 2025 [shift-photocurrent vortices](https://doi.org/10.1038/s41699-025-00552-6) work has a different drive/output; **Introduction/nearest texture precedent**, not proof that the new exact spin response is pre-empted.

**Locally executed group tests**, in [protected-branch Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32):
- Hex \(C_3+m_v+m_z\): real \(\mathbf j_\parallel\to S_-^z\) response **six symmetry-forced zero locations G/A/B/C/M/N**; explicit Fourier witness simple-zero indices \(+1,+1,+1,-1,-1,-1\), sum 0.
- Square \(C_4+m_v+m_z\): **four** forced zero locations; witness indices \(+1,+1,-1,-1\), sum 0. Exact rational fixed-point solver confirms \(|\det(I-C_3)|=3\), \(C_4=2\), layer inversion \(=4\), unions 6 vs 4.
- **Sharper same-system readout contrast:** hex layer-odd \(S_-^z\) has first **linear** slide response \(F_-\sim i\gamma u\) (winding +1 if nonzero), while layer-even \(S_+^z\) has first **quadratic** \(F_+\sim\beta\bar u^2\) (winding −2 if nonzero). A separate symmetry-consistent even Fourier field passes +1/+1/−2 index budget. Square \(C_{2z}\) vs horizontal layer swap forces **all-registry layer-even \(S_+^z\) current response identically zero** (known centrosymmetric-selection consequence).

These are **conditional, standard symmetry/topological deductions**, NOT a new microscopic spin mechanism, completed 80/528 sweep or novelty proof. The actual current-driven electronic coefficients remain uncalculated. The literature-overlap test with the original 2026 LEE SM remains open. Registry [now 28 ideas](data/BD_IDEA_REGISTRY_20261022.json), with BD27 and BD28 still OPEN_HYPOTHESIS; no new DFT/HPC.
