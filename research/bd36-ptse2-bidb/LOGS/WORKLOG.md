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

## BD36 R1.1 corrections (2026-10-23, per user review of Draft PR #99)

1. POSCAR_SD selective dynamics: bilayers x,y frozen / z free (T F F); monolayer
   full relax (T T T). Post-relax verification added (verify_relaxed.py): inversion/
   C3/mirror re-audit + registry drift + slide residual; dry-run PASS (baseline D3d
   retained; slides Ci retained; slide residual ~1e-16 vs declared +-(0.08,0.19)).
2. Wannier projections unified: Pt:d + Se:p only; num_wann = 2*(5*nPt+3*nSe) =
   44 (bilayer) / 22 (monolayer); num_bands=64 >= num_wann.
3. VASP/W90 interface fixed: LWANNIER90=.TRUE. + wannier90.win as FILE; win
   atoms_frac generated per structure from that POSCAR (slide-correct, verified);
   win kpoints NOT hard-coded - make_win_kpoints.py generates post-SCF from the
   actual VASP k-list (ISYM=-1, 9x9x1 -> 81); no embedded WANNIER90_WIN, no
   spin=true. Platform variants: vasp6 (IMR) / vasp544 (para).
4. Spin matrix path: LWRITE_SPN=.TRUE. (VASP-written .spn). para read-only probe:
   wannier90.x == 3.0.0; VaspWannier vasp_ncl binary CONTAINS LWANNIER90/
   LWRITE_SPN/LWRITE_UNK keywords (strings), BUT the HgI2 L1 benchmark has NO
   outputs -> end-to-end .spn NOT validated -> in-place validation job required
   before production. IMR: version UNVERIFIABLE from this sandbox (no pty);
   LWRITE_SPN by VASP>=6.4.3 design, pending user confirmation. read_spn.py added.
5. Total spin (S_i = <w_n|sigma_i|w_m>, .spn) vs layer-projected spin
   (S_layer = P_layer S_z P_layer by Wannier centers) distinguished in PLAN.md
   §5; response via PT-degenerate-subspace projected traces (gauge-safe).
6. Caps fixed: per job 28 cores x <=2 h (<=56 core-h); round-1 total <=500 core-h
   (minimal set 448); committed scripts cleaned of internal absolute paths
   (env placeholders + local run_env.sh, gitignored).
7. static_check_inputs.py: 17 checks/structure + budget; 6/6 structures PASS
   (0 FAIL), exit 0. make_win_kpoints.py + read_spn.py + verify_relaxed.py added.

No production compute submitted. PR #99 updated with per-item PASS/FAIL.
