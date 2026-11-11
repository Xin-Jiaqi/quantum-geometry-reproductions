# BD36 R1 | PtSe2 BiDB layer-Edelstein inputs audit & plan (READ_ONLY)

Branch: `research/bd36-ptse2-bidb-r1-20261023` (isolated; main + other branches untouched).
Draft PR #99 (this branch) builds on PR #98 (inputs kit).

Deliverables in this directory:
- `PLAN.md` - research plan (small-first: 2253 zero-response baseline + (0.08,0.19) slide + 2252 monolayer; deferred: reverse slide, 2255, 2254)
- `STRUCTURE_SYMMETRY_AUDIT.md` - pure-numpy site-permutation symmetry audit (all 6 POSCARs) + selection rules
- `INPUTS/` - frozen, hashed VASP/Wannier input sets (relax / soc / soc+w90 / KPOINTS / POSCAR_SD / wannier90.win / PBS+Slurm drafts, all NOT_RUN) + CONVERGENCE_TABLE.md
- `COST_AND_PERMISSION_GATE.md` - read-only environment inventory + job cost table + approval request
- `scripts/` - audit_structures.py, build_inputs.py (reproducible, no DFT/no SSH/no POTCAR)
- `LOGS/` - run logs

Safety: no qsub/VASP/remote modifications/WSL/Slurm submissions were made. POTCAR and BiDB SQLite are not committed.
