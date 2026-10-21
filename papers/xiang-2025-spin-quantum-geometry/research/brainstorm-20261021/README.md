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
