# Codex entry point | Xiang 2025 spin quantum geometry + separate LEE research

Read this file before making changes. Then read `CODEX_HANDOFF_20261023.md` and the scientific acceptance file `docs/ACCEPTANCE.md`. **Do not treat the adjacent sliding LEE idea as reproduction of the Xiang PRL.**

## Two independent research tracks

**Track A (paper reproduction).** Xiang, Jin & Wang, *Spin Transport Revealed by Spin Quantum Geometry*, PRL 135, 146303 (2025). Frozen paper-model/method implementation in this folder. Scientific status **PARTIAL**, despite tests passing. Models include two-valley Dirac NDSC, reduced Eq. 13, Fig. 1(b) distinct source formulas, SM spin QGT, finite-eta/PT projector method, V31 Laurent and V32 source Eq. 2 checks. Open: Fig 1(a) ≈2× absolute axis, SM Eq 55/56 and Eq 69/70 algebra, metallic driven DC interpretation, source/full DC protocol. Original complete offline V1–V27 archives are **not fully** in this repository. Never silently change frozen reference data or "repair" a printed equation.

**Track B (exploratory physics).** Stacking-dependent layer-odd out-of-plane spin response `j_xy -> S_top^z - S_bottom^z` in identical parallel bilayers, originally in `research/brainstorm-20261021`; BD31 through BD36. Closest existing literature: Zhou et al., *Layer Edelstein Effect* PRL 137, 076401 (2026). An earlier idealized horizontal-mirror family is **not** the identical centrosymmetric 1T-PtSe2 family: carefully distinguish their symmetry and zeros. BD32–BD35 are *model/theory*, not real-material PtSe2 evidence. BD36 is a BiDB structure and compute-environment investigation, not a first-principles spin response. See `CODEX_HANDOFF_20261023.md` and PR #99.

## Local startup (READ_ONLY; no DFT)

```bash
git status --short --branch
python -m pip install -r papers/xiang-2025-spin-quantum-geometry/requirements.txt
cd papers/xiang-2025-spin-quantum-geometry
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

If local shell or repo layout differs, adapt only commands, not scientific acceptance criteria. Python 3.11 recommended. Show exact test count, failures and commit SHA. For the exploratory BD32–BD35 snapshots, use the separate portable archive and checksum listed in `CODEX_HANDOFF_20261023.md` and `research/brainstorm-20261021/tasks/CODEX_LEE_CURRENT_20261023.md`, after the archive is present locally. Do not pretend source archives are committed when they are not.

## Operating constraints

1. Never run an old MnPS3 R1 command from `research/brainstorm-20261021/tasks/CURRENT.md`: R1 is CLOSED, both prior jobs timed out before ionic steps, R2-P0 only READ_ONLY.
2. Never run the old para BD36 Phase-0 script: physical para preflight failed, no .spn-capable binary and no licensed Pt/Se POTCAR, zero jobs submitted; its script was hard-blocked in PR #99. The old one-job authorization applies **only to that now-invalid para job**, not replacement IMR/WSL/QE jobs.
3. No new `sbatch`, `qsub`, DFT production, HPC, remote edits, private POTCAR/credentials uploads or auto-retries without **new explicit bounded user authorization for the exact task**.
4. Before proposing novel LEE mechanisms, label each claim **SYMMETRY THEOREM / MATHEMATICAL WITNESS / UNFITTED TOY / DATABASE STRUCTURE / DFT MATERIAL / NOVELTY OPEN**. Never promote toy numerical zero coordinates to PtSe2 predictions.
5. Preserve existing modules and baselines, work in an isolated branch, add minimal regression tests; report commands and source provenance.
6. Work independently: inspect GitHub PRs directly and write results back to the PR; do not require the user to relay assistant messages.

## Suggested first autonomous actions

- Read `CODEX_HANDOFF_20261023.md`, then inspect the actual GitHub PR #99 `FAST_IDEA_GO_NOGO_20261023.md` and `LOGS/phase0_preflight_report.md`. Only READ_ONLY environment inventory is presently within prior scope.
- Confirm whether BD32–BD35 original source archive has been staged locally. If absent, mark **MISSING_LOCAL_ARCHIVE**, do not invent script contents. Request only one-time archive placement rather than repeated long reports.
- Run available tests, create a compact PASS/FAIL status sheet and stop. Do not initiate another open-ended R1.4/R1.5 review loop.

Last human-approved state: separate branches; no new production compute permission.