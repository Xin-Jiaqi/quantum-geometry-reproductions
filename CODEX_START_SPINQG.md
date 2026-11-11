# Local Codex | spinQG single-entry handoff (2026-10-23)

**Start with this file.** The complete repository work is divided into:
1. **Xiang/Jin/Wang 2025 spin quantum geometry PRL** original paper reproduction, code is already in GitHub, SCIENCE=PARTIAL.
2. **Separate BD31–BD36 layered Edelstein research** (theory, synthetic model, BiDB structures, and blocked para preflight); not the original spinQG paper reproduction.

**Codex instruction:**

> Read `papers/xiang-2025-spin-quantum-geometry/AGENTS.md`, `papers/xiang-2025-spin-quantum-geometry/CODEX_HANDOFF_20261023.md`, and `papers/xiang-2025-spin-quantum-geometry/research/brainstorm-20261021/tasks/CODEX_LEE_CURRENT_20261023.md`. Check the actual Git state and paper scientific acceptance. Execute only local Python paper-model tests and read-only research checks. Cross-check current BD36 para preflight/physics-first plan at PR #99 without running remote jobs. Log evidence to an isolated branch/PR. Do not modify the paper baselines or restart MnPS3/para/VASP jobs. If the BD32–BD35 ZIP is available locally, use `scripts/import_bd32_bd35_local.py` to verify SHA and import to the separate research directory; if missing, mark IMPORT_BLOCKED_MISSING_ARCHIVE rather than guessing. Report only genuinely unsolved scientific and resource decisions.

**One-time local archive**: `spinqg_handoff_sources_bd32_bd35_20261023.zip` (SHA256 `605843ccf73d64d29f05da7181b915429167337a852f2b025f6adfc1c66d1602`). The archive is attached in the ChatGPT conversation, **not in this GitHub branch** until local Codex imports and commits it. There is no automatic way for this GitHub repository to read a ChatGPT sandbox attachment.

**GitHub references:** Xiang spinQG: `papers/xiang-2025-spin-quantum-geometry/` in main; LEE / BiDB + blocked Phase 0: [Draft PR #99](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/99). Old `research/brainstorm-20261021/tasks/CURRENT.md` is **MnPS3 history**, not a current compute authorization.

**Research-integrity status:** no new production DFT, Slurm or PBS jobs were run/authorized by creating this handoff.