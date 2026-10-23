#!/usr/bin/env python3
"""BD36 build_inputs v2: materialize frozen, hashable VASP/Wannier input sets.

READ_ONLY builder. No SSH, no POTCAR handling, no qsub/sbatch, no remote writes.
POTCAR intentionally excluded (private; assembled at execution from the approved
registry, SHA256 recorded locally). All energies/responses NOT_RUN.

R1.1 fixes (per user review of Draft PR #99):
1. POSCAR_SD selective dynamics: bilayers freeze x,y (registry preserved), z free;
   monolayer (2252) relaxes fully (T T T) as the independent control.
2. Wannier projections unified to Pt:d + Se:p; num_wann = 2 x (5*nPt + 3*nSe):
   bilayer 44 spinor, monolayer 22 spinor. No s-projections.
3. VASP/W90 interface: LWANNIER90=.TRUE. + wannier90.win as a FILE (no embedded
   WANNIER90_WIN, no spin=true in win); win atoms_frac generated per structure
   from that structure's POSCAR (slide-correct); win kpoints NOT hard-coded
   (generated post-SCF from the actual VASP k-list by make_win_kpoints.py).
4. Spin matrices: LWRITE_SPN=.TRUE. path (VASP-written .spn), platform variants:
   - IMR/VASP 6.4.x (vasp6 variant) - native; version must be user-confirmed.
   - para/VaspWannier 5.4.4 (vasp544 variant) - keywords verified present in the
     binary, but end-to-end .spn output NOT yet validated -> in-place validation
     job required before production.
5. Total vs layer-projected spin documented in PLAN.md; response uses layer
   projectors P_layer built from Wannier centers, PT-degenerate-subspace traces.
6. Job caps explicit (per job <= 28 cores x 2 h; round-1 total <= 500 core-h);
   committed scripts contain NO internal absolute paths (env placeholders only).
7. static_check_inputs.py verifies projections/k-points/POSCAR-Win consistency/
   tag compatibility/resource budget; verify_relaxed.py re-audits CONTCAR.
"""
import hashlib, json, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "structures"
OUT = ROOT / "INPUTS"

ENCUT = 500
KGRID = "9 9 1"
CAPS = {"cores": 28, "wall_h": 2, "job_core_h": 56, "conv_core_h": 28,
        "phase0_validation_core_h": 56,   # para step-0 in-place interface validation
        "phase1_round1_core_h": 448,       # 3 structures x 2 jobs x 56 + 4 conv x 28
        "total_core_h": 504}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_poscar(path):
    lines = open(path).read().splitlines()
    lat = [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coords = [[float(x) for x in l.split()[:3]] for l in lines[8:8 + sum(counts)]]
    return lines, lat, species, counts, coords


def magmom_line(n):
    return " ".join(["0 0 0.02"] * n)


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
LASPH     = .TRUE.
GGA_COMPAT = .FALSE.
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
LASPH     = .TRUE.
GGA_COMPAT = .FALSE.
"""

INCAR_SOC_W90_VASP6 = """# BD36 SOC + VASP6 native W90 interface (IMR target)
SYSTEM   = {tag} soc w90 interface (vasp6)
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
ISYM     = -1
LWANNIER90 = .TRUE.
LWRITE_SPN  = .TRUE.
LWRITE_UNK  = .FALSE.
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
LASPH     = .TRUE.
GGA_COMPAT = .FALSE.
"""

INCAR_SOC_W90_VASP544 = """# BD36 SOC + VaspWannier(v5.4.4) W90 interface (para variant)
# NOTE: LWANNIER90/LWRITE_SPN keywords verified present in this build's binary,
# but end-to-end .spn output NOT yet validated -> in-place validation job first.
SYSTEM   = {tag} soc w90 interface (vasp544)
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
ISYM     = -1
LWANNIER90 = .TRUE.
LWRITE_SPN  = .TRUE.
LWRITE_UNK  = .FALSE.
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
LASPH     = .TRUE.
GGA_COMPAT = .FALSE.
"""

W90_WIN = """# BD36 {tag} wannier90.win (spinors; .spn via VASP LWRITE_SPN, NOT spin=true)
num_wann = {num_wann}
num_bands = {nbands}
spinors = true
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
Pt:d
Se:p
end projections
begin unit_cell_cart
Ang
{cell}
end unit_cell_cart
begin atoms_frac
{atoms}
end atoms_frac
# begin kpoints: GENERATED POST-SCF from the actual VASP k-list
# (scripts/make_win_kpoints.py IBZKPT|OUTCAR) - do not hand-edit.
"""

KPOINTS = """Gamma-centered {kgrid} (convergence-tested, base mesh)
0
Gamma
{kgrid}
0 0 0
"""

PBS_DRAFT = """#!/bin/bash
# PBS Pro draft (IMR, Japan). NOT_RUN - explicit user approval required.
# Toolchain locations come from the LOCAL env file (run_env.sh, NOT committed).
# No internal absolute paths in this file.
#PBS -N BD36-{tag}
#PBS -q {queue}
#PBS -l select=1:ncpus={cores}:mpiprocs={cores}:mem={mem}gb
#PBS -l walltime={wall_h}:00:00
#PBS -j oe
#PBS -V
set -e
cd "$PBS_O_WORKDIR"
[ -f ./run_env.sh ] || {{ echo "missing run_env.sh (local, uncommitted)"; exit 2; }}
source ./run_env.sh
: "${{VASP_BIN:?}}"
$MPI_RUN -np {cores} "$VASP_BIN" > vasp.stdout
echo "EXIT_OK"
"""

SLURM_DRAFT = """#!/bin/bash
# Slurm draft (para). NOT_RUN - explicit user approval required.
# Toolchain locations come from the LOCAL env file (run_env.sh, NOT committed).
#SBATCH --job-name=BD36-{tag}
#SBATCH --partition=root.bingxing
#SBATCH --nodes=1
#SBATCH --ntasks={cores}
#SBATCH --time={wall_h}:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -eo pipefail
ROOT="${{SLURM_SUBMIT_DIR:?}}"
cd "$ROOT"
[ -f ./run_env.sh ] || {{ echo "missing run_env.sh (local, uncommitted)"; exit 2; }}
source ./run_env.sh
: "${{VASP_BIN:?}}"
$MPI_RUN -np {cores} "$VASP_BIN" > vasp.stdout
echo "EXIT_OK"
"""

RUN_ENV_EXAMPLE = """#!/bin/bash
# LOCAL-ONLY toolchain environment (git-ignored; do NOT commit).
# Fill in per approved platform. No internal paths belong in committed scripts.
export VASP_BIN=          # e.g. path/to/vasp_ncl (SOC build)
export MPI_RUN=mpirun     # e.g. mpirun (Intel MPI) or mpirun.openmpi
# optional for para VaspWannier build:
# export LD_LIBRARY_PATH=...mkl/intel64:...mpi/lib:...
"""

STRUCTURES = {
    "2252_monolayer": ("POSCAR_BiDB_2252_PtSe2", "round1", "monolayer"),
    "2253_bilayer_baseline": ("POSCAR_BiDB_2253_PtSe2", "round1", "bilayer"),
    "2253_generic_p008_p019": ("POSCAR_PtSe2_BiDB2253_generic_p008_p019_rigid", "round1", "bilayer"),
    "2253_generic_m008_m019": ("POSCAR_PtSe2_BiDB2253_generic_m008_m019_rigid", "round1_deferred_job", "bilayer"),
    "2255_bilayer": ("POSCAR_BiDB_2255_PtSe2", "deferred", "bilayer"),
    "2254_bilayer": ("POSCAR_BiDB_2254_PtSe2", "deferred", "bilayer"),
}


def sd_flags(n, kind):
    # VASP selective-dynamics semantics: T = coordinate FREE, F = coordinate
    # FIXED. Bilayers: x,y FIXED (F F) + z FREE (T) -> "F F T", preserving the
    # in-plane registry exactly. Monolayer: full relaxation "T T T".
    if kind == "monolayer":
        return ["T T T"] * n          # full relaxation (independent control)
    return ["F F T"] * n              # registry-preserving: x,y FIXED, z FREE


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    manifest = {}
    for tag, (posname, phase, kind) in STRUCTURES.items():
        src = STRUCT / posname
        lines, lat, species, counts, coords = read_poscar(src)
        n = sum(counts)
        d = OUT / tag
        d.mkdir()
        shutil.copy(src, d / "POSCAR")
        flags = sd_flags(n, kind)
        with open(d / "POSCAR_SD", "w") as f:
            f.write("\n".join(lines[:7]) + "\nSelective dynamics\n" + lines[7] + "\n")
            for ln, fl in zip(lines[8:8 + n], flags):
                f.write(ln.rstrip() + "  " + fl + "\n")
        npt = sum(c for s, c in zip(species, counts) if s == "Pt")
        nse = sum(c for s, c in zip(species, counts) if s == "Se")
        num_wann = 2 * (5 * npt + 3 * nse)
        cell = "\n".join("  " + " ".join(f"{v: .12f}" for v in row) for row in lat)
        sp_list = [s for s, cnt in zip(species, counts) for _ in range(cnt)]
        atoms = "\n".join(
            f"{sp} " + " ".join(f"{v: .12f}" for v in c)
            for sp, c in zip(sp_list, coords))
        win = W90_WIN.format(tag=tag, num_wann=num_wann, nbands=64, mp_grid=KGRID,
                             cell=cell, atoms=atoms)
        (d / "wannier90.win").write_text(win)
        (d / "INCAR_relax").write_text(INCAR_RELAX.format(tag=tag, encut=ENCUT))
        magmom = magmom_line(n)
        (d / "INCAR_soc").write_text(INCAR_SOC.format(tag=tag, encut=ENCUT, magmom=magmom))
        (d / "INCAR_soc_w90_vasp6").write_text(
            INCAR_SOC_W90_VASP6.format(tag=tag, encut=ENCUT, magmom=magmom, nbands=64))
        (d / "INCAR_soc_w90_vasp544").write_text(
            INCAR_SOC_W90_VASP544.format(tag=tag, encut=ENCUT, magmom=magmom, nbands=64))
        (d / "KPOINTS").write_text(KPOINTS.format(kgrid=KGRID))
        (d / "PBS_draft.pbs").write_text(PBS_DRAFT.format(
            tag=tag, cores=CAPS["cores"], mem=128, wall_h=CAPS["wall_h"], queue="P_030"))
        (d / "sbatch_draft.sh").write_text(SLURM_DRAFT.format(
            tag=tag, cores=CAPS["cores"], wall_h=CAPS["wall_h"]))
        (d / "run_env.example.sh").write_text(RUN_ENV_EXAMPLE)
        (d / "README").write_text(
            f"BD36 {tag} ({phase}; {kind}). Inputs frozen, NOT_RUN.\n"
            f"POSCAR_SD selective-dynamics semantics (VASP): T = coordinate FREE, F = coordinate FIXED.\n"
            f"POSCAR_SD: {'monolayer full relax T T T (all coordinates FREE, no fixed components)' if kind == 'monolayer' else 'bilayer registry-preserving F F T (x,y FIXED - in-plane registry exactly preserved; z FREE - interlayer relaxation only)'}.\n"
            f"wannier90.win: num_wann={num_wann} (Pt:d + Se:p x spinor), atoms_frac from THIS POSCAR "
            f"(slide-correct); kpoints block GENERATED post-SCF by make_win_kpoints.py.\n"
            f"Interface: LWANNIER90 + LWRITE_SPN (vasp6 or vasp544 variant); no spin=true, no embedded WANNIER90_WIN.\n"
            f"POTCAR absent (private; assembled at execution, SHA256 recorded locally).\n"
            f"Caps: {CAPS['cores']} cores x {CAPS['wall_h']} h per job; phased: phase0 validation {CAPS['phase0_validation_core_h']}, phase1 round-1 {CAPS['phase1_round1_core_h']}, total {CAPS['total_core_h']} core-h.\n")
        manifest[tag] = {
            "poscar": posname, "poscar_sha256": sha(src), "phase": phase, "kind": kind,
            "num_wann": num_wann, "n_atoms": n,
            "files": sorted(p.name for p in d.iterdir()),
        }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    (OUT / "caps.json").write_text(json.dumps(CAPS, indent=2, sort_keys=True))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
