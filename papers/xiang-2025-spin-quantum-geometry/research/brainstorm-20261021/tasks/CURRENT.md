# CURRENT TASK | B+D real-material research

**Last updated:** 2026-10-22.  
**Active science-review task:** **R2-P0 READ-ONLY APPROVED / NOT CONFIRMED STARTED BY LOCAL HARNESS**. **Compute execution:** NO ACTIVE JOBS. R1 remains CLOSED; no R2 VASP/WSL/para/Japan submission authorized.

## Most recent task: MnPS3 L2 polar-pair R1 — CLOSED

- **Task:** BD-MNPS3-L2-P1-R1
- **Source:** bilayer MnPS3 published P1+(L2) ↔ P1-(L2) path 1 (Sun et al., Nat. Commun. 17, 3101, 2026).
- **Decision:** accepted final R1 outcome, and **declined additional submission**.
- **Preflight:** source-derived structural/magnetic mapping PASS; L2, Cm, published Eq.13 geometric relation checked, but starting structures were **not verified author-relaxed structures**.
- **Real execution:** two authorized paired para Slurm ionic-relaxation attempts; both reached two-hour scheduler TIMEOUT and were terminated. Scheduler exit 0:0 is **not** numerical success.
- **Scientific acceptance:** RUN_FINISHED=true (terminal job lifecycle), NUMERIC_CONVERGED=false, FE_VERIFIED=false; 0 ionic steps, no validated total-energy difference, polarization, spin transport or torque data.
- **Public closeout with limitations:** [MnPS3 R1 scientific status and user decision](../docs/BD_MNPS3_R1_CLOSEOUT_20261022.md).
- **Private source of truth:** authenticated \`Xin-Jiaqi/wannauto-research-archive\` R1 \`result.json\` / \`results_summary.md\` in \`research/20261021/bd-mnps3-l2-p1-r1/\`. Private full fingerprints, queue evidence and input files remain there or in compute storage.
- **Historical R1 instruction (DO NOT REEXECUTE):** [BD-MNPS3-L2-P1-EXECUTE-R1.md](BD-MNPS3-L2-P1-EXECUTE-R1.md).

## Directive for a local DeepSeek harness

**Do not execute the old R1 command just because you can read this public CURRENT pointer.** R1 is already closed. Do **not** re-submit interrupted ionic jobs, switch hosts, raise parallelism to 64 ranks, start Wannier or run an R2 campaign without a *new task ID* and the user's explicit additional scientific/resource approval.

If a user asks for *status*, read the private dated acceptance manifests and actual live scheduler when accessible; report the final R1 closeout above rather than stale \`READY\` or \`RUNNING\` text from a historical instruction. If a fresh R2 command is issued, use the archived R1 SCF oscillation/divergence evidence and frozen inputs to formulate **a controlled low-cost stabilization test**, not an unbounded rerun, and perform preflight before any submission.

The former single-link \`CURRENT.md\` URL still works, but it is a **state pointer**, not a background runner. It will only authorize execution again if a new task is explicitly recorded and separately approved.

## Sources

- [Real-material R0 scientific preflight](../docs/BD_MNPS3_REAL_MATERIAL_PREFLIGHT_20261021.md)
- [Physics and local harness handoff contract](../docs/BD_RESEARCH_TO_COMPUTE_HANDOFF_20261021.md)
- [Independent Gate 5 toy Ward response](../docs/BD_GATE5_TORQUE_WARD_20261021.md)
- [Scientific literature review](../docs/PEER_WORK_AUDIT_20261021.md)

Do not upload private POTCAR, system paths, large complete matrix data, private unpublished structures or raw complete job logs to this public repository.


## Approved zero-compute next task: MnPS3 R2-P0

**Versioned, actionable task:** [BD-MNPS3-L2-R2P0-READONLY.md](BD-MNPS3-L2-R2P0-READONLY.md). The user has explicitly approved this **read-only analysis** in conversation. It asks the local DeepSeek harness to autonomously inspect **existing private R1 four-Mn site/moment mapping and frozen SCF failure evidence**, run [four-site operator tests](../scripts/mnps3_l2_site_operator_gate.py), distinguish the author's **nonrelativistic [E||MzMy]** from any physical **SOC magnetic symmetry**, define a source-correct Néel-mode generalized force, and return a PRIVATE audit with a sanitized status and **proposed** one-variable SCF stabilization pilots.

**New 2025–2026 papers, reviewed:** [Guo vdW SOT, Haque MnPS3 THz, Ilyas FePS3 THz, Kang octupolar AFM tunnel torque](../docs/BD_PEER_UPDATE_20261022_SOT_THz_OCTUPOLE.md), now integrated into the cumulative peer audit.

**Authorizations:** Scientific read-only review YES. **New calculation/spending/Slurm/WSL/Japan/para jobs NO.** R1's two TIMEOUTs are final; any future R2 compute would need a **separate task and user authorization**. A GitHub document does **not** itself start the remote local harness: no R2-P0 execution receipt or private site-map acceptance has been received as of this public entry.

### Local harness single-link instruction

> 读取本 CURRENT.md 中已批准的 R2-P0 READ-ONLY 任务。自动检查现有 R1 私有 MnPS3 L2 +/− 原始位点、磁序和日志；运行已归档的四站点算符代数测试和 SOC 对称性审核，生成私有验收及公开脱敏摘要。不要让我搬运材料；不得调用 VASP、sbatch 或重新启动任何 R1 作业。若必须依赖缺失的结构证据，明确标记 BLOCKED，而不是猜测。

**Research source of truth:** private R1 \`result.json\` for historical job completion and nonconvergence; R2-P0 output **not yet received**. No claim of material magnetic response or new effect.


## Scientific priority note (2026-10-22)

[**High-return physics ideas, peer work and falsification/kill gates**](../docs/BD_IDEA_SEARCH_DECISION_20261022.md): Prioritize an experimentally distinguishable **ferroelectric path-dependent AFM normal-mode drive** in a source-valid MnPS3 L2 site basis. The read-only active R2-P0 should report \`FE state/path × real Neel mode × symmetry-allowed E field × closest SOT precedent\`; it may assess a **second path at source/group-algebra level only**. No new DFT, Wannier, NEB, job-array or Slurm task is permitted by this note. If all relations reduce to existing tensor/normal-mode selection rules, **close this idea as NO-GO** without a speculative materials calculation.


## Physics focus pivot: group-wide symmetry, not more MnPS3 DFT (2026-10-22)

The user requested a **high-value general layer-group/representation theory research pivot**. Source: [full physics/novelty/implementation decision memo](../docs/BD_GROUP_THEORY_80_528_PIVOT_20261022.md); group-toolkit protected-main [Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32). Current proved algebraic negative control: exact UNITARY full magnetic partner domains **cannot generate a new independent mode-driving channel** merely by FE domain reorientation. This is useful standard equivariance, not independently novel.

**Next scientific G1:** find one source-compatible **magnetic-site mode representation** and a FE switching operation **not** conjugating the FULL SOC magnetic Hamiltonian, alongside an exactly conjugate negative-control pair. Compare stabilizer \`Hom_H(E,mode)\` dimensions and physical mode mapping before any 80/528 table or DFT. If every apparent bright/dark case is just rotated drive axes, call NO-GO. See new 2026 spin-layer-group / magnon / FE prior-art constraints in the cumulative peer report.

The previously approved **MnPS3 R2-P0 remains a separate READ-ONLY private evidence task**, not confirmed received; it can supply actual ordered sites but is **not a dependency for starting abstract group analysis**. **R1 CLOSED, no R2 compute/spending permissions or Slurm/WSL/Japan re-submission.**


## 2026-10-22 science-review follow-up; unchanged compute permissions

[**Full original-paper / brainstorm / peer-work audit**](../docs/BD_CROSS_PROJECT_DEEP_AUDIT_20261022.md) supersedes overly optimistic interpretation of Gate1–6 and PR #32. The established **antiunitary finite-frequency Onsager reciprocity involves exchanged operators**, unlike the previously tested **unitary static covariance**. This mathematical correction has been tested in a standalone two-/four-band Lehmann [group-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32).

**Research priority next:** G2 representative site-/magnetic-mode representations and honest unitary/antiunitary/nonconjugate domain gate; optional moire local-response texture is a separate research hypothesis. **No new MnPS3, WSL, para or Japan jobs permitted.** MnPS3 R1 remains CLOSED and source-real R2-P0 remains **APPROVED READ-ONLY / execution not confirmed**.


## Group theory G2 result and G3 scientific criterion (2026-10-22)

**[G2 p3 layer-group Wyckoff/sublattice gate](../docs/BD_G2_P3_SITE_ORBIT_FROBENIUS_20261022.md)** is complete as **unitary Γ point group algebra**, archived in [group-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) with a site-stabilizer Frobenius reciprocity crosscheck. It derives a **real** contrast of site-local longitudinal \`E→S_z\` allowed-space dimension 4 vs 0 for six-site generic vs special positions, with **no macroscopic signal in either**, and an independent transverse spin 8 vs 12 dimension. This is standard, already researched local hidden response symmetry, **NOT a known switchable FE pair** or new material measurement. Other rejection gates: \`rankK≤2\` for uniform in-plane E and real \`C3: KᵀK=λI\`.

**G3 search priority:** physically legitimate FE magnetic-connection graph with site orbits/mode projectors, not an arbitrary Wyckoff-occupancy comparison; distinguish unitary full-domain covariance, antiunitary **operator-swapped Onsager**, and genuine SOC-magnetically inequivalent switches. Test closest 2018–2026 multipoles/hidden-spin competitors *before* starting all-LG/MLG enumeration. \`R1\` remains CLOSED and the separate MnPS3 \`R2-P0\` is **READ-ONLY**, with **ZERO** authorized Slurm/WSL/Japan VASP jobs.


## 2026-10-22 new idea-falsification archive and research boundary

**Scientific framework:** [8-stage Idea Kill Gates, 20 preserved claim verdicts](../docs/BD_IDEA_ELIMINATION_GATES_20261022.md) and [JSON registry](../data/BD_IDEA_REGISTRY_20261022.json). [Group-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) contains a conservative static idea preflight and example fixtures; not a registered production API or novelty oracle, exact Python execution **not verified here**.

**G3 next hypothesis, not yet an ACTIVE compute task:** [real sliding FE domain-wall / magnetic-site mode force, specific closest peer checks](../docs/BD_G3_NEXT_HYPOTHESIS_20261022.md). Distinguish **real FE wall dynamics** from endpoint Seitz connection, ordinary local SOT/multipoles and electric bilayer torque. No real source-backed magnetic FE wall or material Kubo has been accepted; status is OPEN_HYPOTHESIS / BLOCKED_EVIDENCE, not reported discovery.

**Compute control unchanged:** R1 CLOSED; R2-P0 remains independently approved READ-ONLY; no Slurm, WSL, para, Japan or new VASP task is authorized by this addition. Reading this note does not dispatch a local DeepSeek task.


## 2026-10-22 G3 cost decision — no wall DFT; method-only registry gradient

The user raised DFT cost/uncertainty of a real domain. [**G3 ROI and scientific preemption report**](../docs/BD_G3_DOMAIN_DFT_COST_AND_REGISTRY_GRADIENT_20261022.md): full magnetic FE-wall DFT/torque is **DEFERRED**; even the two primitive MnPS3 R1 SCFs reached walltime without convergence. Existing FE wall Edelstein, magnetic sliding bilayer spin/orbital Edelstein, phonon-driven Edelstein and actual sliding-mode coherent dynamics block broad novelty. [24-claim idea registry](../data/BD_IDEA_REGISTRY_20261022.json) updated.

**New zero-DFT method prototype** in [group toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32): registry-gradient \`Hom_H(U⊗E,O)\` vs clamped \`Hom_H(E,O)\` with explicit site/relative-layer representation warning. **Scientific next gate, not a compute authorization:** verify \`U\` for two source-backed BSF cases, one with and one without layer exchange; test physical response novelty before sample SCF. The old MnPS3 R2-P0 is independently READ_ONLY and R1 remains CLOSED. No Slurm, VASP, WSL, Japan or additional billing authorized.


## Physics-first seed & claim-relative references (2026-10-22)

**[Core critique, direct-competitor vs Introduction rules, and one worked symmetry seed](../docs/BD_CORE_PROBLEM_REFERENCE_ROLES_AND_SMALL_SEED_20261022.md).** User's research correction: build one small real response first; do not automatically escalate from group selection rules to 528 groups, domain-wall DFT or a complete paper outline. \`NO_GO_CLAIM\` is limited to the exact claim, NOT a topic-wide prohibition. For the specifically chosen current-induced layer-odd Sz registry derivative, [2026 Layer Edelstein PRL](https://doi.org/10.1103/qhvs-pfkq) is the **closest direct research framework**, whereas magnetic local torque, WTe2 sliding phonon, and FE-wall magnetism are **context** unless we choose their actual same observable/protocol.

**New completed non-DFT test:** a three-generator ideal bilayer (C3z, vertical mx, horizontal layer-swapping mz) admits precisely one \`relative-slide u × in-plane current j → layer-odd Sz\` FIRST-order invariant, while forbidding its zero-slide term; 4×4 constraint SVD checked. This is **not** a new physical mechanism; gamma coefficient unknown. Test scripts archived in [group-theory PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32).

**Next scientific step (not a compute task):** inspect the original LEE source SM equation-by-equation for same \`d(chi_S,j)/du\`, then calculate a single controlled \`chi(u)\` curve in a symmetry-compatible SOC/transport model. No DFT/HPC/Slurm authorized. MnPS3 R1 CLOSED, separate R2-P0 still READ-ONLY, not auto-dispatched.


## 2026-10-22 BD25 exact model test and next evidence gate

[Source-model audit and dimensionless Boltzmann baseline](../docs/BD25_LEE_MODEL_AUDIT_20261022.md): the published LEE Eq. (3) has no \(S_z\) texture, hence \(\chi_{S_-^z,j_y}=0\) identically in its band-diagonal Boltzmann regime; scalar registry-dependent \(t_\perp\) or \(\Delta\) cannot test the earlier proposed slide derivative. This narrows the problem rather than invalidating the symmetry result. BD25 now has status BLOCKED_EVIDENCE: it needs a source-backed spinful model with both \(S_z\) matrix elements and microscopic \(\partial_u H\). Do not add a free spinful hopping merely to force nonzero \(\gamma\). Original LEE PRL is the closest framework; its SM still needs equation-level overlap review. No DFT/Slurm/WSL/Japan task.


**Model search update (22 Oct):** the [BD25 model audit](../docs/BD25_LEE_MODEL_AUDIT_20261022.md) found a possible source-backed route: 11-orbital MoS\(_2\) Slater–Koster bilayer TB with registry-dependent S–S/Mo–S/Mo–Mo hopping (PRB 102, 081103(R), 2020) plus a separately sourced SOC extension. Not implemented yet. Next gate is only the \(u=0,\pm\delta\) real-space hopping/Seitz covariance test, not Kubo/DFT. PRB 114, 125129 (2026) and CrI\(_3\) 2024 pre-empt generic sliding-tunable Edelstein novelty but are strong-neighbor/Introduction context for the narrower nonmagnetic \(S_-^z\) derivative. Registry now has 26 claims.


**Latest low-cost test:** [AA MoS\(_2\) registry geometry gate](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/aa_mos2_registry_geometry_gate.py) passed locally: 450 interlayer pair/translation checks under \(m_z\) and \(\mathbf u\to-\mathbf u\), max distance mismatch \(1.776\times10^{-15}\) Å; S–S interface spacing 3.6763 Å matches the source SK model. **Geometry only**. Next gate is the actual orbital Slater–Koster matrix covariance including orbital parity/rotation; still no SOC/Kubo or production compute.


## 2026-10-22: abstraction-level correction and general-theory hypothesis

User's key correction: a specific sliding-ferroelectricity example does not automatically pre-empt a general bilayer-stacking theory; a specific ferroic phase does not automatically pre-empt a general symmetry principle. Prior art must be matched to the **scope of the exact proposition**, not merely the phenomenon. See [scope-aware literature policy and candidate theory](../docs/BD_THEORY_LEVELS_AND_STACKING_RESPONSE_FIELD_20261022.md). All 27 registry entries now carry working scope-level and contribution-type labels.

**BD27 (OPEN_HYPOTHESIS, S2–S3):** symmetry-covariant layer-resolved response field over the stacking-registry torus. For a layer-odd \(S_z\) driven by in-plane current, a layer-swapping horizontal-mirror Seitz operation gives \(\chi(-\delta)=-\chi(\delta)\); candidate consequence: zeros at symmetry-fixed registries \(2\delta\in\Lambda\) and forbidden even Taylor orders. This is ordinary covariance unless it yields a nontrivial global result not already in prior literature. Next: verify actual Seitz operations, compare exact scope with LEE SM and stacking-configuration theories. No DFT/HPC authorized.


**Prior-art check for BD27:** [Ji et al., PRL 130, 146801 (2023)](https://doi.org/10.1103/PhysRevLett.130.146801) explicitly generalize earlier sliding-ferroelectric examples into an 80-layer-group BSF theory; [Pang & He, PRL 135, 116402 (2025)](https://doi.org/10.1103/trhd-kxm1) generalize FQFE and conventional ferroelectricity under a generalized Neumann principle. These are the right examples of different abstraction levels. They also mean BD27 must not claim a new general Neumann principle or generic stacking-space symmetry classification. Its only plausible opening is a nontrivial theorem for a specified driven layer-resolved response field that is not already implied by BSF, LEE, or the generalized-Neumann framework.


## 2026-10-22 BD27/BD28: tested general response-field seed, zero DFT

**[Complete mathematics/physics gate](../docs/BD27_REGISTRY_TORUS_ZERO_INDEX_GATE_20261022.md)** now derives and numerically tests a cross-family operator-dependent response field on the stacking-registry torus:
- Hexagonal \(C_3+m_v+m_z\): **six** symmetry-forced zeros G/A/B/C/M/N for \(\mathbf j_\parallel\to S_-^z\), vs only **three** \(C_3\)-forced zeros for \(S_+^z\). Near G, lowest allowed \(S_-^z\) term is **linear**, index +1 if nonzero; \(S_+^z\) is **quadratic**, index −2 if nonzero.
- Square \(C_4+m_v+m_z\): **four** symmetry-forced layer-odd zeros, while the **layer-even** \(S_+^z\) current response is globally zero for all registry positions via \(C_{2z}\) and horizontal layer-swap symmetry.
- Independent Fourier witnesses and exact rational \((I-R)\delta\in\Lambda\) fixed-point code passed locally on the [group-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32); vortex charge budgets match torus Euler number 0.

**Crucial evidence status:** textbook covariance / Poincaré–Hopf and specific high-symmetry BSF points are already known; this is a **candidate for an S2 operator-resolved global selection classification**, not a newly claimed theorem, real Kubo result or experimental prediction. LEE SM comparison remains open. Same-observable direct competitor determination must happen at the precise theorem/response level, not by the existence of sliding-FE examples or different shift-current vortices. Registry now [28 claim entries](../data/BD_IDEA_REGISTRY_20261022.json).

**Next high-value action:** compare the 2026 LEE SM on the exact \(S_\pm^z\), current protocol and global registry fields, then derive a second source-compatible real layer/site response; only after a microscopic nonzero coefficient is shown consider a wider 80 LG theory. R1 CLOSED, R2-P0 read-only, **no Slurm/DFT/WSL/Japan authorization**.
