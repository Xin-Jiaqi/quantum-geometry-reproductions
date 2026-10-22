# AgF₂ G1 · SCF Audit (pending — job queued, not started)

**Queue state as reported by the local Harness, timestamp pending independent reconciliation**: job `490620.cmspbs1` has been submitted
exactly once (no auto-retry, no scope expansion) and is **still queued** on
`CP_001` (`job_state = Q`; PBS comment "Not Running: Insufficient amount of
resource: Qlist"). No OUTCAR/vasprun.xml exists yet, so **no SCF result is
claimed**. This audit file records the *plan* and will be updated with real
values when the job completes.

**Provenance warning:** the original receipt used `2026-10-21T04:05:00Z` as its as-of time, but its companion input record reports submission at `2026-10-21T04:14:00Z`. The chronology is internally inconsistent. Both original timestamps remain preserved as reported values in JSON and must be checked against actual PBS `qstat -f` / accounting timestamps. No live queue state has been independently observed by the GitHub reviewer, and no VASP result is yet available.

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

---

# Update 2026-10-21 (later): local WSL diagnostic run — REAL SCF evidence (user-requested)

**Not** the approved PBS execution (`490620.cmspbs1` is still queued; that
receipt's facts above remain untouched). The user explicitly requested running
the identical single point on their local WSL machine (private VASP 6.6.0
build, OpenMPI 4.1.6, 8 ranks). Input bytes are identical to the PBS
submission (all four SHA-256 match, including the locked POTCAR
`8c0da262…`). Full facts: `G1_LOCAL_WSL_DIAGNOSTIC.json`.

| item | value (local run, vasp.6.6.0 06Mar2026 complex) |
| --- | --- |
| EDIFF=1E-6 reached | **yes** — OUTCAR: `aborting loop because EDIFF is reached` |
| Davidson iterations | 27 |
| TOTEN (eV / 12-atom cell) | **−35.17579601** (single point, unrelaxed geometry) |
| E-fermi (eV) | −5.6235 |
| net total magnetization (μB) | 0.0000 (OSZICAR `mag=-0.0000`) |
| local Ag moments (μB, POSCAR order) | **−0.567, +0.567, +0.567, −0.567** (d-channel ±0.575) |
| F moments | \|m\| ≤ 0.001 (noise) |
| magnetic pattern | Type-I interlayer Néel OPPOSITE seed `[-1,+1,+1,-1]` **preserved**; compensated AFM-like state (signs relative, no SOC easy axis) |
| preliminary gap (12×8×1, occ≥0.5) | ≈ **1.49 eV** (HOMO −5.9015 / LUMO −4.4076) — NOT k-converged, diagnostic only |
| CPU / elapsed | 1193.8 s / 1194.8 s (~20 min), 8 ranks, peak ~3.2 GB |
| error/interruption | none; clean timing/accounting block, exit 0 |

Operational note: the first local launch attempt died silently during setup
when the controlling ssh session closed (no SCF output); relaunched under a
persistent session and completed cleanly. Full outputs (OUTCAR, vasprun.xml,
vaspout.h5, OSZICAR, EIGENVAL, DOSCAR, PROCAR) retained privately on the WSL
host.

Guardrails: single point at the UNRELAXED C2DB-derived geometry; 12×8×1 and
ISMEAR=0/SIGMA=0.05 are screening choices; the ~1.49 eV gap and ±0.57 μB
moments are **preliminary diagnostics**, not a reproduction claim. The PBS
job (vasp.6.4.3) will supply the version-matched counterpart when it runs;
this local run does not change the queued-only PBS receipt above.
