# BD36 R1 worklog (2026-10-23)

- Read DEEPSEEK_TASK.md + PR #98 inputs (6 POSCARs + manifest + README).
- Built pure-numpy site-permutation symmetry audit; results match BiDB manifest
  (layer groups 72 p-3m1 / 69 p3m1; slides retain inversion, break C3+mirrors;
  permutation errors ~1e-8 A; slide shift measured = declared +-(0.08,0.19),
  layer gap preserved). See STRUCTURE_SYMMETRY_AUDIT.md + LOGS/audit_run.log.
- Environment gate (read-only): para = verified SOC+W90 chain (vasp_ncl,
  wannier90.x, Pt/Se POTCAR_REC, HgI2 L1 template); IMR versions unverifiable
  from sandbox (pty); WSL has Se POTCAR, Pt missing, VASP6+postw90 present.
- Materialized INPUTS for 6 structures (3 job types each + KPOINTS + POSCAR_SD
  + wannier90.win drafts + PBS/Slurm drafts), all NOT_RUN, hashed in manifest.
- Wrote PLAN.md / STRUCTURE_SYMMETRY_AUDIT.md / COST_AND_PERMISSION_GATE.md /
  INPUTS/CONVERGENCE_TABLE.md.
- No production compute submitted (permission gate).
