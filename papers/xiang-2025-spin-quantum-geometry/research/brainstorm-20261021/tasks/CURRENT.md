# CURRENT TASK | B+D real-material research

**Last updated:** 2026-10-22.  
**Current execution status:** **NO ACTIVE COMPUTE TASK; R1 CLOSED BY USER**. No R2 jobs authorized, no auto-retry, no active submission requested here.

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


## Read-only follow-on prepared (not active computing)

[**R2-P0 site- and mode-resolved magnetic operator audit — DRAFT, zero compute**](BD-MNPS3-R2P0-SITE-OPERATORS-DRAFT.md) is now available to the local DeepSeek harness **if the user explicitly asks it to perform the review**. It requires reading existing private R1 source/site mapping, distinguishing honeycomb Mn Néel from toy layer-odd spin, testing the actual SOC magnetic tensor and planning one-variable-at-a-time SCF stabilization **without any new Slurm/WSL/Japan compute job**.

**R1 remains CLOSED**, with no new authorization to submit or restart jobs. Do not reinterpret this note as a launch signal.
