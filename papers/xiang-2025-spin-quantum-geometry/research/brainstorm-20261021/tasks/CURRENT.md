# CURRENT TASK | B+D material computation handoff

**Updated:** 2026-10-21 (research task specification only; **NOT an actual remote queue status**).  
**Canonical task ID:** BD-MNPS3-L2-P1-R1  
**Current mission:** [Execute R1: MnPS₃ L₂ published FE pair, verify sources and submit at most two bounded WSL/Japan baseline jobs](BD-MNPS3-L2-P1-EXECUTE-R1.md).

## For the local DeepSeek / Codex harness

**Treat this GitHub file as the stable task entry.** Always fetch the latest version from its actual URL. Follow the linked, versioned R1 instruction *in order*, including all STOP gates and 2-job resource cap. Do not infer from the status text on GitHub that jobs have been submitted or that remote hosts are accessible. Use the host's live queue, verified input/commit hashes and local private acceptance manifests as truth. Only proceed to DFT once original or source-justified structures and precise L₂ magnetic mappings are accepted.

**Do not require the user to paste separate R0 reports, VASP INCAR, command outputs or logs.** Fetch the linked public reports yourself and inspect the *authorized* local repositories and existing host state. Source paper is [Sun et al., Nat. Commun. 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2). The correct L₂ path-1 source operation is **author Eq.13 [E || M_z M_y]**, NOT L₁ Eq.11 or previous toy-model m_z T.

Status in this public task is only **AVAILABLE_FOR_AGENT_INTAKE**. The actual live execution state must be established by the local agent from its own private manifest and platform job IDs. This file does not run a GitHub Action or remotely contact WSL or Japan on its own.

## Reusable activation message

> 请直接读取 GitHub 中的 CURRENT 任务入口以及其中引用的 R1 完整执行指令，自己获取文献结构、检查本地 DeepSeek/WSL/日本超算既有配置，先完成零成本结构和磁序预检；证据充分且权限、预算都满足时自动提交**最多两项**低成本配对 DFT 基准任务，追踪并进行数值验收。不要让我搬运文件或日志；一旦遇到缺失原始结构、未授权资源或无法确认的物理设置，停在明确的 BLOCKED 状态，只问一个合并问题。不要修改原论文复现，也不要将私有原始数据提交到公开 GitHub。

## Research provenance

- Prior source and detailed science: [published-structure R0 review](../docs/BD_MNPS3_REAL_MATERIAL_PREFLIGHT_20261021.md).
- Existing local-execution contract: [GPT–DeepSeek–WSL/HPC scientific handoff](../docs/BD_RESEARCH_TO_COMPUTE_HANDOFF_20261021.md).
- Toy-only result: [Gate5 spin-transfer/torque Ward experiment](../docs/BD_GATE5_TORQUE_WARD_20261021.md).
- Peer-work landscape: [literature audit](../docs/PEER_WORK_AUDIT_20261021.md).

No secret token, scientific unpublished structure or live queue information is placed here.
