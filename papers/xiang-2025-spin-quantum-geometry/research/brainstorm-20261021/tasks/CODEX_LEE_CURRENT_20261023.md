# Codex current task: spinQG context + PtSe2 layer-Edelstein exploration

**Date:** 2026-10-23. **Task type:** local research handoff / read-only auditing. **Status:** READY_FOR_LOCAL_INTAKE; material computation NOT AUTHORIZED.

**One canonical instruction:** Read `papers/xiang-2025-spin-quantum-geometry/AGENTS.md` and `papers/xiang-2025-spin-quantum-geometry/CODEX_HANDOFF_20261023.md` in the `handoff/spinqg-codex-20261023` branch. The old `tasks/CURRENT.md` is the MnPS3 R1-closed / R2 read-only historical pointer, **not this LEE task**.

## Facts already confirmed

1. Paper reproduction Xiang 2025 Spin Quantum Geometry: curated runnable code in GitHub `main`; software tests ≠ full article acceptance, paper science PARTIAL.
2. BD31 idealized registry symmetry and exact Seitz tests already on GitHub or group-theory PR #32.
3. BD32–BD35 analyses currently exist as independently downloaded 43KB portable ZIP, **not yet committed into GitHub**. They have a local source-import helper and hashes on this handoff branch. Import before claiming complete transfer; do not ask user to manually paste source code.
4. BD36 structural/compute preflight in separate PR #99; para actual preflight FAIL (no Pt/Se POTCAR; no .spn writer); **zero new job ID / zero core hours**. Old para Phase-0 script hard blocked. Phase 1 not approved.
5. New physics-first staged goal is G0 platform readiness, then one general-slide SOC layer-projected spin doublet trace, then coarse DC current ratio F=Q D^{-1}, only then detailed Wannier. G1/G2/G3 are **NOT RUN**.

## Actions local Codex can take immediately without additional resource approval

- Perform `git status`, inspect existing paper module and verified scientific acceptance, run normal Python paper-model unit tests.
- If user has placed the research ZIP locally, run `python papers/xiang-2025-spin-quantum-geometry/scripts/import_bd32_bd35_local.py --archive /absolute/path/spinqg_handoff_sources_bd32_bd35_20261023.zip`. The importer checks the archived SHA256, prevents overwrites, and writes a file manifest. Run BD32 site and torus scripts locally; log return codes. More expensive BD33–BD35 runs can be separately scoped to avoid runaway resource use.
- Read both `research/bd36-ptse2-bidb/LOGS/phase0_preflight_report.md` and `FAST_IDEA_GO_NOGO_20261023.md` from PR #99 without merging it, mark G0 platform INCONCLUSIVE until environment read-only audit locates appropriate SOC binary and properly licensed Pt/Se potentials. Do NOT attempt a new production job.
- Open a draft PR from an isolated local working branch containing imported research code, source hash manifest, offline Python test logs, and precise remaining physics blockers. Link PR #99 rather than duplicating old blocked compute inputs.

## Stop and report rather than improvise

- Missing source ZIP: report **IMPORT_BLOCKED_MISSING_ARCHIVE** and the exact expected filename/hash.
- Missing licensed Pt/Se potentials: report **G0_BLOCKED_INPUT**, do not borrow others' potentials or make computational substitutes.
- No trustworthy layer-spin operator: report **G1_INCONCLUSIVE**, do not interpret random PT Kramers band spin projections as a physical observable.
- Any new compute estimate: request explicit authorization of platform, job number, runtime, memory, structures, outputs. Earlier 56 core-hour *para* authorization is not transferable.
- Preserve original source-paper reproduction and old MnPS3 closed tasks.

## Minimal acceptance record

```text
GIT_COMMIT=
PAPER_MODEL_TESTS=NOT_RUN/PASS/FAIL (counts/log path)
BD32_SOURCE_IMPORT=NOT_PRESENT/IMPORTED + sha256
BD32_SEITZ_TEST=NOT_RUN/PASS/FAIL
BD32_TORUS_WITNESS_TEST=NOT_RUN/PASS/FAIL
BD35_SYNTHETIC_LOCKING_CAVEAT=ACKNOWLEDGED
BD36_PARA=BLOCKED_PREFLIGHT (no job)
G0_ENV=NOT_RUN/INCONCLUSIVE/PASS
G1_LAYER_SPIN=NOT_RUN
G2_EDELSTEIN_RATIO=NOT_RUN
G3_FULL_WANNIER=NOT_RUN
MATERIAL_VALIDATION=NO
NEXT_ONE_DECISION=
```

The next action is **source migration and read-only audit**, not another theoretical Fourier toy and not VASP production. Record all status in GitHub, reducing user involvement to actual choices and permissions.