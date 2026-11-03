# DeepSeek Harness task: one SOC-on AgF₂ B1 static screening point

**Created 2026-10-22. Status: PROPOSED_NOT_EXECUTED. This file itself is not an execution receipt or proof of compute authorization.** The researcher should explicitly send the short authorization message at the end before the Harness starts.

## Goal

Answer one cheap, high-information question: does the screened AgF₂ B1 Type-I candidate support a self-consistent SOC-on insulating magnetic state compatible with the intended x-axis Néel seed, or does SOC materially change its local moments, gap, or magnetic symmetry?

This is a **single-point SOC screen only**. It does not compute IGMC, magnetic anisotropy energy, a magnetic ground state, structural stability, Berry polarization, B2, NEB, Wannier interpolation, or a response amplitude.

## Source of truth

Repository: https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38  
Branch: research/agf2-b1-wsl-relax-20261021  
Source geometry: idea_screen_1/POSCAR.screen_last_complete  
Prior run receipt: idea_screen_1/RUN_RECEIPT.json  
Prior report: idea_screen_1/IDEA_SCREENING_VERDICT.md  
Tensor correction / research motivation:
- idea_screen_1/B1_IGMC_TENSOR_ERRATUM_20261022.md
- idea_screen_1/B1_BSF_IGMC_NOVELTY_STORY_GATE_20261022.md

The source structure is a C2DB-derived screened candidate, not the authors' relaxed POSCAR. Use the **exact submitted step-11 coordinates** and verify their file hash before copying. Do not silently re-relax, alter the cell, change atom order, or substitute a different geometry.

The previous non-SOC result: 11 ionic steps; max force 0.0782 eV/Å at screening threshold; local Ag moments about [-,+,+,-]; occupation-defined gap 1.3875 eV on 12×8×1. It is not a proven local minimum or magnetic ground state.

## Mandatory preflight; if any gate fails, mark BLOCKED and do not run

1. Check the existing local VASP installation has an executable **vasp_ncl** compatible with the installed VASP 6.6.0 + OpenMPI environment. **LSORBIT=.TRUE.** automatically enables noncollinear spinors and requires vasp_ncl. **Do not build/recompile VASP, install another version, switch to vasp_std, or submit PBS.**
2. Check memory, disk, MPI launch and the prior-run directory. Prior screening used 8 ranks and about 7.8 GB WSL-visible RAM; SOC can cost more. Use at most 8 MPI ranks and do not start if the current machine cannot safely support the job.
3. Create a **new isolated directory** under the private local run root. Copy only POSCAR, INCAR, KPOINTS and the licensed private POTCAR. Do not copy WAVECAR/CHGCAR; the task intentionally starts from scratch. Never overwrite prior G1/B1 files.
4. Confirm POSCAR atom order is Ag4 F8 and the concatenated Ag/F POTCAR hash matches the locked private hash **8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff**. Never upload POTCAR or full private output files.

## Exact single-point settings

Use the same 12×8×1 Gamma KPOINTS and screened B1 geometry. No ionic steps.

~~~text
SYSTEM = AgF2_B1_TypeI_SOC_static_screen_1
GGA = PE
PREC = Accurate
ENCUT = 500
EDIFF = 1E-6
NELM = 160
NELMIN = 4
ALGO = Normal

LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = 2 -1
LDAUU = 4 0
LDAUJ = 0 0
LMAXMIX = 4
LASPH = .TRUE.
IVDW = 11

LSORBIT = .TRUE.
LNONCOLLINEAR = .TRUE.
SAXIS = 0 0 1
MAGMOM = -1 0 0  1 0 0  1 0 0  -1 0 0  24*0.0
ISYM = -1
GGA_COMPAT = .FALSE.
LORBIT = 11
LORBMOM = .TRUE.

ISMEAR = 0
SIGMA = 0.05
ISTART = 0
ICHARG = 2
LREAL = .FALSE.
NCORE = 4
NSW = 0
IBRION = -1
ISIF = 2
LWAVE = .FALSE.
LCHARG = .FALSE.
~~~

### Input rationale and pitfalls

- There must be **no ISPIN=2** line. For noncollinear calculations VASP ignores ISPIN, and VASP 6.5+ can error if ISPIN=2 is combined with LNONCOLLINEAR.
- MAGMOM has three components per ion: four Ag moments seed the previously selected Type-I [-,+,+,-] pattern along Cartesian x; all F moments start at zero. Default SAXIS=(0,0,1) makes spinor axes coincide with Cartesian axes in the documented VASP convention. Record the exact SAXIS and MAGMOM.
- MAGMOM is an initialization, **not a constraint**. Do not add I_CONSTRAINED_M. If the self-consistent spin directions rotate, report the actual vectors; do not automatically restart along y or z.
- ISYM=-1 avoids symmetry symmetrization in this first SOC screen. Classify the converged magnetic point group afterward from the final structure and spin vectors; do not simply stamp m'_x onto the result.
- GGA_COMPAT=.FALSE. and LASPH=.TRUE. are recommended for noncollinear GGA precision.
- LORBMOM requests site/orbital angular-momentum diagnostics. These do not by themselves give orbital IGMC.
- Keep EDIFF=1e-6, Ueff(Ag-4d)=4 eV, D3 zero-damping (IVDW=11), ENCUT and k mesh fixed for this first screen. Do not attempt to optimize k convergence or MAE in this task.

## Resource limit and execution policy

- One local WSL static SCF only, maximum 8 MPI ranks, maximum 4 hours wall-clock (32 rank-hours upper cap).
- No Japan HPC/PBS, no second launch, no retry, no parameter backoff, no ionic relaxation, no second spin orientation, no Type-II/strained branch, no B2, no Wannier, no IGMC response calculation.
- If vasp_ncl is missing, memory is unsafe, the launcher fails, or the job cannot be monitored safely, report **BLOCKED**; do not fix the environment by compiling or expanding permissions.
- If the wall-clock cap is reached, report the actual partial status and preserve the final complete electronic iteration if possible. Do not silently exceed the cap. Do not pre-create STOPCAR.
- Do not claim a converged result unless OUTCAR explicitly reaches EDIFF and a complete final electronic energy/force/eigenvalue/magnetization record exists.

## Required output in the same PR #38

Add under a new idea_screen_1/soc_static_1/ subdirectory:

1. RUN_RECEIPT.json: status (PASS_SCREEN, PARTIAL, FAIL, or BLOCKED), actual VASP binary/version, ranks, walltime, source and input hashes, stop reason, SCF convergence, and whether the task stayed inside scope.
2. SOC_STATIC_VERDICT.md: concise scientific interpretation, no more than 2 pages.
3. Sanitized INCAR.soc, KPOINTS.soc, and a copy of the exact screened POSCAR as POSCAR.soc_input; no POTCAR or full private raw logs.
4. MOMENTS.csv: site-indexed spin moment vector and orbital moment vector if reliably present; distinguish spinor basis from Cartesian basis and record SAXIS. Missing values must remain null/blank, never invented.
5. BANDS_GAP.json or equivalent: occupation-defined gap method, HOMO/LUMO, k mesh and clear warning that the gap is only a screening value.
6. A concise SCF diagnostic / selected log markers, plus SHA256 hashes of private raw outputs in the receipt. Raw OUTCAR/vasprun.xml/POTCAR remain private.
7. Comment in PR Conversation with status, final SCF energy, gap/metallic diagnosis, local spin vectors, whether the conditional m'_x hypothesis survives, and exactly one next recommendation.

## Interpretation rules

- PASS_SCREEN: SCF converges and the state remains physically interpretable (local compensated AFM moments and a clear screening gap). This is **not** a proof of the easy axis, MAE, magnetic ground state, stable local minimum, ferroelectric switching, or nonzero IGMC.
- PARTIAL: SCF is incomplete or the gap / magnetic state cannot be established. No automatic retry.
- FAIL: actual converged result shows the prepared magnetic seed collapses/reconstructs materially or the candidate is metallic under this setup. This is a result for this seed/approximation, not proof the material can never support IGMC.
- BLOCKED: prerequisites fail or safe single-run execution cannot be guaranteed.

A successful static point only answers whether this SOC-included electronic state merits a next-stage **response-theory feasibility check**. It does not authorize any further calculation.

## Official VASP guidance

- LSORBIT and required executable / symmetry caution: https://vasp.at/wiki/LSORBIT
- Noncollinear MAGMOM basis and spinor setup: https://vasp.at/wiki/MAGMOM
- Noncollinear executable and ISPIN restriction: https://vasp.at/wiki/LNONCOLLINEAR

## One-message authorization prompt

After reviewing this task, the researcher can send the Harness:

“Please read PR #38 task file B1_SOC_STATIC_SCREENING_TASK_20261022.md and execute exactly the one SOC-on static WSL screening point it specifies, if and only if all preflight gates pass. I authorize only that single run within the 8-rank / 4-hour cap. Do not submit PBS or perform retries, relaxation, MAE orientations, B2, Wannier, or IGMC. Commit sanitized results directly to the same PR and stop.”
