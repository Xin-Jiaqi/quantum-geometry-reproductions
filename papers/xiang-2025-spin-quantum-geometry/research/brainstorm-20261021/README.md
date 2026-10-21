# Independent A/B/C/D brainstorming (2026-10-21)

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
