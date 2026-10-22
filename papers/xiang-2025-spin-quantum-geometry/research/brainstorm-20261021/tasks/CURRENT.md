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
