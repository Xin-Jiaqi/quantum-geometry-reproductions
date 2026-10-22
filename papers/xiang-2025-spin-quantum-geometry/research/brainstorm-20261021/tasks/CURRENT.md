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
