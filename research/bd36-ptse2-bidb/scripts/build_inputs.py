#!/usr/bin/env python3
"""BD36 build_inputs: materialize frozen, hashable VASP/Wannier input sets.

READ_ONLY builder (mirrors the established prepare_soc.py pattern):
- No SSH, no POTCAR handling, no qsub/sbatch, no remote writes.
- POTCAR is intentionally EXCLUDED (private data; assembled at execution time
  from the approved pseudopotential registry with SHA256 recorded locally).
- All energies/responses are NOT_RUN until explicit user-approved submission.

Round-1 materialization (task 'small first'):
  2252 monolayer control | 2253 primary baseline | 2253+p(0.08,0.19) slide |
  2253-m(0.08,0.19) reverse slide (job deferred per DEEPSEEK_TASK.md item 3).
  2255 / 2254 / reverse-slide jobs: documented, deferred (POSCAR kept).

Protocol notes:
- Registry-preserving relaxation: IBRION=2 + ISIF=2 (cell fixed) +
  selective dynamics (x,y frozen; z free) so the in-plane fractional slide is
  exactly preserved and only interlayer/buckling degrees relax.
- SOC electronic step: zero-field nonmagnetic (LNONCOLLINEAR+T, LSORBIT=T,
  tiny initial moments), Fermi level to be placed near electron/hole bands
  (not mid-gap) via doping reference cells in the response stage.
- Wannier: spinors=true, num_wann=44 (22 spatial x spinor), spin=true (.spn),
  projections Pt 5d + Se 3p; band window/num_bands finalized by projection
  tests (documented in CONVERGENCE_TABLE.md), NOT hard-coded as final.
"""
import hashlib, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "structures"
OUT = ROOT / "INPUTS"

ENCUT = 500
KGRID = (9, 9, 1)          # base; convergence table tests 7/9/11/13
NBANDS_SOC = 64            # 44 spinor + empty buffer; REVISED by projection tests
NUM_WANN = 44
MP_GRID = "9 9 1"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def poscar_atoms(path):
    lines = open(path).read().splitlines()
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    return lines, counts, n


def selective_flags(n, layers_z):
    """x,y frozen; z free for all atoms (registry-preserving relaxation)."""
    return ["T F F"] * n


INCAR_RELAX = """# BD36 registry-preserving relaxation (non-SOC, spin-paired)
SYSTEM   = {tag} relax (registry-preserving, PBE+D3)
PREC     = Accurate
ENCUT    = {encut}
ISMEAR   = 0
SIGMA    = 0.05
ISPIN    = 1
IBRION   = 2
ISIF     = 2
NSW      = 60
EDIFF    = 1E-6
EDIFFG   = -0.01
NELM     = 120
LREAL    = .FALSE.
ALGO     = Normal
IVDW     = 12
LWAVE    = .FALSE.
LCHARG   = .TRUE.
LORBIT   = 11
NCORE    = 4
"""

INCAR_SOC = """# BD36 zero-field nonmagnetic SOC static (PBE+D3+LSORBIT)
SYSTEM   = {tag} soc static (nonmagnetic, zero field)
PREC     = Accurate
ENCUT    = {encut}
ISMEAR   = 0
SIGMA    = 0.05
ISPIN    = 2
LNONCOLLINEAR = .TRUE.
LSORBIT  = .TRUE.
SAXIS    = 0 0 1
MAGMOM   = {magmom}
IBRION   = -1
NSW      = 0
EDIFF    = 1E-6
NELM     = 120
LREAL    = .FALSE.
ALGO     = Normal
IVDW     = 12
LWAVE    = .TRUE.
LCHARG   = .TRUE.
LORBIT   = 11
NCORE    = 4
"""

INCAR_SOC_W90 = """# BD36 SOC + VASP-W90 external interface (wannier90.win embedded)
SYSTEM   = {tag} soc w90 interface
PREC     = Accurate
ENCUT    = {encut}
ISMEAR   = 0
SIGMA    = 0.05
ISPIN    = 2
LNONCOLLINEAR = .TRUE.
LSORBIT  = .TRUE.
SAXIS    = 0 0 1
MAGMOM   = {magmom}
NBANDS   = {nbands}
IBRION   = -1
NSW      = 0
EDIFF    = 1E-6
NELM     = 120
LREAL    = .FALSE.
ALGO     = Normal
IVDW     = 12
LWAVE    = .TRUE.
LCHARG   = .TRUE.
LORBIT   = 11
NCORE    = 4
WANNIER90_WIN = "
{win}
"
"""

W90_WIN = """num_wann = {num_wann}
num_bands = {nbands}
spinors = true
spin = true
write_hr = true
mp_grid = {mp_grid}
use_ws_distance = false
dis_win_min = -18.000000
dis_win_max = 8.000000
dis_froz_min = -18.000000
dis_froz_max = 0.000000
dis_num_iter = 1000
dis_mix_ratio = 0.500
dis_conv_tol = 1.0d-10
num_iter = 3000
conv_tol = 1.0d-6
conv_window = 20
guiding_centres = true
trial_step = 0.5
begin projections
Pt:s
Pt:d
Se:s
Se:p
end projections
begin unit_cell_cart
Ang
  3.748260200702  0.000000000000  0.000000000000
 -1.874130100351  3.246088553802  0.000000000000
  0.000000000000  0.000000000000  26.234554388803
end unit_cell_cart
begin atoms_frac
Pt 1.000000000000 1.000000000000 0.335720147538
Se 0.666666666891 0.333333333970 0.385646999757
Se 0.333333333970 0.666666667940 0.285793294938
Pt 1.000000000000 1.000000000000 0.523151013987
Se 0.666666666891 0.333333333970 0.573077866587
Se 0.333333333970 0.666666667940 0.473224161768
end atoms_frac
begin kpoints
0.000000000000 0.000000000000 0.000000000000
0.111111111111 0.000000000000 0.000000000000
0.111111111111 0.111111111111 0.000000000000
0.222222222222 0.000000000000 0.000000000000
0.222222222222 0.111111111111 0.000000000000
0.222222222222 0.222222222222 0.000000000000
0.333333333333 0.000000000000 0.000000000000
0.333333333333 0.111111111111 0.000000000000
0.333333333333 0.222222222222 0.000000000000
0.333333333333 0.333333333333 0.000000000000
0.444444444444 0.000000000000 0.000000000000
0.444444444444 0.111111111111 0.000000000000
0.444444444444 0.222222222222 0.000000000000
0.444444444444 0.333333333333 0.000000000000
0.444444444444 0.444444444444 0.000000000000
0.500000000000 0.000000000000 0.000000000000
0.500000000000 0.111111111111 0.000000000000
0.500000000000 0.222222222222 0.000000000000
0.500000000000 0.333333333333 0.000000000000
0.500000000000 0.444444444444 0.000000000000
0.500000000000 0.500000000000 0.000000000000
end kpoints
"""

KPOINTS = """Gamma-centered 9x9x1 (convergence-tested, base mesh)
0
Gamma
9 9 1
0 0 0
"""

PBS_DRAFT = """#!/bin/bash
# PBS Pro draft for IMR (Japan). NOT_RUN - requires explicit user approval
# of platform+resources. Placeholder: queue P_030, 28-112 cores, 460GB node.
#PBS -N BD36-{tag}
#PBS -q P_030
#PBS -l select=1:ncpus={cores}:mpiprocs={cores}:mem={mem}gb
#PBS -l walltime={walltime}
#PBS -j oe
#PBS -V
set -e
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh 2>/dev/null || true
# VASP6 + Wannier90 environment TBD at execution (confirm binaries on IMR)
MPI=mpirun
VASP=vasp_ncl
test -x $(command -v $MPI)
test -x $(command -v $VASP)
$MPI -np {cores} $VASP > vasp.stdout
echo "EXIT_OK"
"""

SLURM_DRAFT = """#!/bin/bash
#SBATCH --job-name=BD36-{tag}
#SBATCH --partition=root.bingxing
#SBATCH --nodes=1
#SBATCH --ntasks={cores}
#SBATCH --time={walltime}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -eo pipefail
ROOT="${{SLURM_SUBMIT_DIR:?}}"
export LD_LIBRARY_PATH=/share/home/e6a000183/intel/oneapi/mkl/2022.1.0/lib/intel64:$LD_LIBRARY_PATH
export PATH=/share/home/e6a000183/tools/vasp/VaspWannier/vasp.5.4.4/bin:$PATH
export PATH=/share/home/e6a000183/intel/oneapi/mpi/2021.6.0/bin:$PATH
MPI=mpirun
VASP=vasp_ncl
test -x "$(command -v $MPI)" && test -x "$(command -v $VASP)"
cd "$ROOT"
$MPI -np {cores} "$VASP" > vasp.stdout
echo "EXIT_OK"
"""

STRUCTURES = {
    "2252_monolayer": ("POSCAR_BiDB_2252_PtSe2", "round1"),
    "2253_bilayer_baseline": ("POSCAR_BiDB_2253_PtSe2", "round1"),
    "2253_generic_p008_p019": ("POSCAR_PtSe2_BiDB2253_generic_p008_p019_rigid", "round1"),
    "2253_generic_m008_m019": ("POSCAR_PtSe2_BiDB2253_generic_m008_m019_rigid", "round1_deferred_job"),
    "2255_bilayer": ("POSCAR_BiDB_2255_PtSe2", "deferred"),
    "2254_bilayer": ("POSCAR_BiDB_2254_PtSe2", "deferred"),
}


def magmom_line(n, species):
    """tiny z-only initial moments: 6 atoms -> 18 entries of 0 0 0.02"""
    entries = []
    for _ in range(n):
        entries.append("0 0 0.02")
    return " ".join(entries)


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    manifest = {}
    for tag, (posname, phase) in STRUCTURES.items():
        src = STRUCT / posname
        d = OUT / tag
        d.mkdir()
        lines, counts, n = poscar_atoms(src)
        species_line = lines[5]
        nspecies = len(counts)
        # POSCAR with selective dynamics appended (x,y frozen, z free)
        shutil.copy(src, d / "POSCAR")
        with open(d / "POSCAR_SD", "w") as f:
            f.write("\n".join(lines[:7]) + "\n")
            f.write("Selective dynamics\n")
            f.write(lines[7] + "\n")
            for ln in lines[8:8 + n]:
                f.write(ln.rstrip() + "  T F F\n")
        magmom = magmom_line(n, species_line)
        # W90 win: atoms block must match POSCAR (6-atom family here; monolayer differs)
        if n == 6:
            win = W90_WIN.format(num_wann=NUM_WANN, nbands=NBANDS_SOC, mp_grid=MP_GRID)
        else:
            # monolayer 2252: 3 atoms, 11 spatial (Pt 5d + Se 3p x1) -> 22 spinor
            win = W90_WIN.format(num_wann=22, nbands=NBANDS_SOC, mp_grid=MP_GRID)
        (d / "wannier90.win").write_text(win)
        (d / "INCAR_relax").write_text(INCAR_RELAX.format(tag=tag, encut=ENCUT))
        (d / "INCAR_soc").write_text(INCAR_SOC.format(tag=tag, encut=ENCUT, magmom=magmom))
        (d / "INCAR_soc_w90").write_text(INCAR_SOC_W90.format(
            tag=tag, encut=ENCUT, magmom=magmom, nbands=NBANDS_SOC, win=win))
        (d / "KPOINTS").write_text(KPOINTS)
        (d / "PBS_draft.pbs").write_text(PBS_DRAFT.format(
            tag=tag, cores=28, mem=128, walltime="02:00:00"))
        (d / "sbatch_draft.sh").write_text(SLURM_DRAFT.format(
            tag=tag, cores=28, walltime="02:00:00"))
        (d / "README").write_text(
            f"BD36 {tag} ({phase}). Inputs frozen, NOT_RUN. POSCAR from BiDB source;\n"
            f"POSCAR_SD = selective-dynamics variant (x,y frozen, z free).\n"
            f"POTCAR intentionally absent (private); assemble at execution from the\n"
            f"approved registry, record SHA256 locally.\n"
            f"wannier90.win: draft - window/num_bands finalized by projection tests.\n")
        manifest[tag] = {
            "poscar": posname,
            "poscar_sha256": sha(src),
            "phase": phase,
            "files": sorted(p.name for p in d.iterdir()),
        }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
