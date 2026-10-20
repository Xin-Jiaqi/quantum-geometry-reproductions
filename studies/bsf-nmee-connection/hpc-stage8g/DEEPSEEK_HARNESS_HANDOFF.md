# Stage 8g | DeepSeek harness autonomous Japan PBS handoff

## Operational requirement

You are the user's already available LOCAL DeepSeek harness with shell, Git, SSH and filesystem abilities. **Execute the work yourself; do not give the user a list of terminal commands to paste.** The user handles scientific choices and access approval, not file transfer, scheduler scripting or log collation.

**Read TASK_CONTRACT.json in this same directory before doing any work.** All existing NMEE paper reproduction and previously frozen reference data are out of scope.

## Known starting point

- Code repository: https://github.com/Xin-Jiaqi/quantum-geometry-reproductions
- Research branch: research/bsf-nmee-stage7-sw-gate-20261020
- Directory: studies/bsf-nmee-connection/hpc-stage8g
- Source geometry: user-uploaded original BiDB row 2685, a WSe2 six-atom z-scan structure already committed in the independent study at material-screening/source-structures/WSe2_BiDB_row2685_3R_polar_zscan.vasp.
- Target protocol: initial **non-SOC, PBE+D3(BJ) 3R-WSe2 AB→BA geometry test**, fixed lattice, freeze each atom's xy and relax z only; this is NOT a rigorous NEB or an NMEE run.
- Previously used Japan resource: **PBS**, with historical DP_002 and CP_001 queues and compute root /work/sunq/xjq. These facts need fresh on-host confirmation; never assume these queues remain currently available. WSL Ubuntu is a known local workspace, NOT necessarily a second HPC scheduler.
- Existing private Wannier workflow: https://github.com/Xin-Jiaqi/wannier-workflow may contain previously validated Japan PBS/oneAPI load scripts. Read/adapt only as needed; do not modify its running jobs or source.

## Automatic workflow: do all executable steps

1. Clone/fetch the named research branch and record HEAD. Work in a newly created directory away from existing AutoWann/NLO work. Verify source geometry and stage8g Python files. Check no scientific source-paper reproduction changes are introduced.
2. From the studies/bsf-nmee-connection directory, run prepare.py in hpc-stage8g with --out pointing to a fresh unique directory OUTSIDE the Git working tree. This produces 15 reproducible path input sets and a SHA-locked manifest. Verify only three are marked pilot: path_A_s0.00, path_A_s0.25 and path_A_s0.50. No POTCAR is included.
3. Inspect WSL and Japan access READ ONLY using existing configured SSH profiles/keys. On Japan check live PBS qsub/qstat, queues, project quota, suitable MPI, licensed VASP, authorized PBE PAW W/Se POTCAR and bootstrap environment. Do not print, store, transmit or commit private key/password/token contents.
4. Prefer Japan PBS if reachable. Transfer the generated job directories, manifest and onsite_pbs.py through the existing authorized SSH/SFTP path. Use an independent run directory, preferably under an existing user-owned runs folder. Avoid using WSL as a job scheduler unless it really has a verified VASP/scheduler stack.
5. Create a **local-only** site.json on Japan. Required keys: scheduler=pbs, queue, ncpus, mem_gb, walltime, vasp_binary (real absolute path), mpi_command (mpirun/mpiexec), license_confirmed and validated. Optional bootstrap_script must be an existing absolute path. Mark validated=true ONLY after probing real executable, permissions, scheduler and queue. Do not push local site.json into a public GitHub repository.
6. Place the locally licensed PBE PAW W and Se potentials in the three pilot directories in the correct species order. Do not transmit POTCAR to GitHub/chat. Verify their hashes and titles match, check ENCUT=520 eV ≥1.3×max(ENMAX), that INCAR/KPOINTS/POTCAR are consistent for all pilot jobs, no structural xy shifts occur during relaxation, and the slab's dipole/damping settings are understood.
7. Execute **dry run first** using onsite_pbs.py --jobs WORKDIR --site LOCAL_SITE_JSON. If status is DRY_RUN_VALIDATED, run the same command with --submit. Maximum **three** first-wave PBS jobs; do not submit the other 12, change inputs silently or rerun already submitted IDs. Save returned IDs in workdir/pbs_ids.json. If environment or authentication cannot pass, report one precise blocker with saved resume state rather than fabricating values.
8. After submission, monitor qstat and preserve remote work/output paths. If your current agent process cannot run persistently, stop safely with actual job IDs and provide the next self-resumable check; never claim finished runs. Once authentic OUTCAR/CONTCAR appear, run onsite_pbs.py --jobs WORKDIR --audit, then independently inspect electronic/ionic convergence, final six-atom forces (max |F_z| preferably <0.01 eV/Å), fixed lattice, true selective constraints, final symmetry and convergence sensitivity. The built-in audit checks final six-atom force tables, z-force threshold, cell invariance and frozen lateral coordinates. Even CONVERGED_REQUIRES_SCIENTIFIC_REVIEW is not automatic scientific acceptance; explicitly inspect electronic convergence and input convergence.
9. Summarize in Chinese: verified site, Git SHA, three PBS IDs and states, remote locations, source input and potential hashes (not potential data), actual total energies if truly run and converged, final max forces, whether symmetry/geometry maintained, missing convergence checks and the next scientific gate. Do not claim CI-NEB barrier, switchable ferroelectricity, SOC insulating full BZ, static NMEE or driven DC.
10. Safe code/provenance corrections may be committed only to the independent research study on the existing Draft PR branch. **Do not** commit output data, POTCAR, credentials or private site settings to a public repo without an explicit data disclosure decision. Leave the source-paper reproduction frozen.

## Abort states (not user busywork)

ABORT and mark BLOCKED if Japanese SSH is unavailable, licensed VASP/POTCAR cannot be proven locally, the configured PBS queue is not verified, an input hash differs, the run directory is not clean, ENMAX requires a convergence adjustment, or any other scientific gate is not satisfied. Describe the exact missing fact and how the harness would resume; ask the user only for access approval or a true scientific decision that cannot be recovered programmatically.

## Physics safeguards

The AB source is an authentic BiDB **z-scan/frozen-monolayer geometry**, not a fully relaxed bilayer. Stage8f established a glide at ideal s=0.5 which forbids the TOTAL in-plane gapped static NMEE, while layer-staggered static response is symmetry-allowed. This is a testable symmetry hypothesis, NOT calculated WSe2 NMEE. The experimental driven-DC Edelstein coefficient in Jia et al. is a distinct limit and is out of scope.

## Final output contract

Use the fields in TASK_CONTRACT.json; report **BLOCKED / SUBMITTED / RUNNING / CONVERGED_REVIEWED**. All empty scientific outcomes remain null. Create next-step requests only when genuinely blocked; the user is not a manual file-moving assistant.
