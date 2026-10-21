# AgF₂ G1 · SCF Audit (pending — job queued, not started)

**As of 2026-10-21 ~04:05 UTC**: job `490620.cmspbs1` has been submitted
exactly once (no auto-retry, no scope expansion) and is **still queued** on
`CP_001` (`job_state = Q`; PBS comment "Not Running: Insufficient amount of
resource: Qlist"). No OUTCAR/vasprun.xml exists yet, so **no SCF result is
claimed**. This audit file records the *plan* and will be updated with real
values when the job completes.

## What will be audited from the real output (private copies stay on the cluster)
| item | source | status now |
| --- | --- | --- |
| VASP exact version | first lines of OUTCAR | pending |
| EDIFF=1E-6 reached | OUTCAR "reached required accuracy" | pending |
| electronic iterations count | OUTCAR step lines | pending |
| TOTEN (eV / 12-atom cell) | OUTCAR "free energy TOTEN" (last step) | pending |
| E-fermi | OUTCAR | pending |
| total magnetization (μB) | OUTCAR "total magnetization" | pending |
| local Ag moments (μB) | OUTCAR "magnetization (x)" block, LORBIT=11 | pending |
| preliminary occupancy/gap | vasprun.xml eigenvalues + occupations; Gaussian smearing caveat | pending |
| CPU/elapsed time | OUTCAR timing block | pending |
| error/interruption | OUTCAR/queue history; Exit_status from qstat -x | pending |
| PBS job history | qstat -x / PBS accounting on the cluster | pending |

## Pre-submission verification already completed (real, recorded)
- POSCAR (GitHub main, frozen) sha256 `2d323c3f…60a` — matched; Ag4F8, Ag/F order, 20 Å vacuum, not author-relaxed.
- INCAR sha256 `8ca091bb…` = `one_shot_g1/INCAR.smoke`; single point `IBRION=-1 NSW=0`; `ISPIN=2`; `MAGMOM=-1 1 1 -1 8*0` (Type-I Néel OPPOSITE seed); `LDAU U_eff(Ag-4d)=4 eV`; `IVDW=11` (D3 zero); `ISMEAR=0 SIGMA=0.05`; `ENCUT=500`; `EDIFF=1E-6`.
- KPOINTS sha256 `4ab9ad06…` = 12×8×1 Gamma (screening, not k-converged).
- POTCAR: assembled privately on the cluster, element order Ag→F; concatenated sha256 `8c0da262…fff` — matched; first TITEL `PAW_PBE Ag 02Apr2005`. Never uploaded.
- VASP binary: `vasp.6.4.3_wannier90_cpu/bin/vasp_std` (the only accessible CPU build for this account; 6.6.1/6.6.0/6.5.1 group-permission denied); ldd after site loader: 0 unresolved; Intel MPI 2021.10.0; launch `mpirun -np 32 -ppn 32 -hostfile $PBS_NODEFILE`.
- Request: 32 ncpus / 120 gb / 02:00:00 on CP_001 — within queue maxima (56/230gb/72h) and the user-approved one-shot budget (64 requested core-hours).

## Boundary statement
`AUTHORIZED_ONE_SHOT` ≠ `physics_approved_by_user=true` for the full Stage-6
scientific scope. `SCF_CONVERGED` (when it happens) is not a reproduction of a
band gap, magnetic ground state, or double well; those require independent
review of moments, occupations, gap and numerical stability.
