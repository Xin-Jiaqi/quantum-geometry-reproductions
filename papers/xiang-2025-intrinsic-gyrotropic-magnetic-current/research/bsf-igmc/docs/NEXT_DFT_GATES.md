# Next-stage gates for AgF₂ BSF × IGMC

The project has moved beyond a purely blocked preflight: one C2DB-derived B1 candidate has passed a low-cost **non-SOC** screening run. The real-material IGMC response remains uncomputed.

## Gate 0 — already passed for screening only

- One local WSL VASP 6.6.0 run, 8 ranks × 1.62 h, 11/11 electronic steps converged.
- RMM-DIIS stopped at max force 0.0782 eV/Å under the deliberately loose screening criterion.
- Type-I local Ag moments remained approximately [-,+,+,-], net moment approximately zero; sampled gap 1.3875 eV; nearest interlayer F–F 3.002 Å.
- These data do not establish a local minimum, SOC easy axis, magnetic ground state, ferroelectric double well, polarization or IGMC.

See PR #38 and the committed screening receipt. Do not rerun this gate.

## Gate 1 — proposed, not executed: one SOC-on B1 static screen

Use the screened B1 coordinates, noncollinear VASP executable, declared x-axis Type-I MAGMOM seed, no ionic relaxation, and a single bounded local WSL run. Verify SCF, site-resolved spin/orbital moments, occupation gap and the actual magnetic symmetry. Full preflight, INCAR and output contract are in [B1_SOC_STATIC_SCREENING_TASK_20261022.md](../agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/B1_SOC_STATIC_SCREENING_TASK_20261022.md). This requires explicit authorization; it is not a task execution record.

## Gate 2 — zero-DFT response-symmetry check

- Distinguish the conduction IGMC tensor C from the displacement IGMC tensor D.
- The earlier static alpha=P/B matrix is D-like only; do not use it as the full IGMC tensor.
- Reconcile the source-level transformation convention with Xiang et al. PRL 134, 116301 and Supplemental Material; retain the project's alternative Eq. (14) convention as an explicit unresolved comparison until audited.
- Use the actual SOC magnetic group and the same laboratory axes for both stacking domains.
- Type-I unitary m_z and Type-II antiunitary m'_z predict different D-channel parity; Type-II is a source-reported strained branch, not an assumed zero-strain state.

See the [tensor erratum](../agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/B1_IGMC_TENSOR_ERRATUM_20261022.md).

## Gate 3 — only if Gates 1–2 pass: quantitative response feasibility

- Source-ground B1/B2 partner structures; distinguish literature-reported structural data from the local C2DB-derived candidate.
- Build SOC Wannier model with same-gauge position and physical spin/Zeeman operators; keep spin-Zeeman and orbital magnetic coupling separate.
- Compute low-cost coarse \(\sigma^D_{xx}(\omega)\) first, with gap/frequency/k-mesh and operator-gauge checks. Compute conduction \(\sigma^C\) only in a physically relevant chemical-potential regime.
- Stop if the response is symmetry forbidden, numerically negligible, or not gauge/convergence stable. No automatic escalation.

## Gate 4 — later, separately authorized

Only after a response signal survives should we invest in high-accuracy relaxation, B1/B2 energy landscape and barrier, Berry-phase polarization, Type-II strained branch, k-mesh/broadening convergence, and publication-grade spin/orbital decomposition.

## Current status

**B1 non-SOC screening: PASS. SOC magnetic state: NOT TESTED. Channel-resolved real-material IGMC: NOT COMPUTED. B1/B2 local minima and material-specific switching barrier in this independent workflow: NOT VERIFIED.** Published AgF₂ BSF results remain literature evidence, not newly reproduced results.

The proposed next task is the single SOC static screen in PR #38; no VASP, PBS, Wannier or response run was started by this document.
