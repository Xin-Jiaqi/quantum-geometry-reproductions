#!/usr/bin/env python3
"""Stage8g GitHub-native, reproducible 3R-WSe2 initial structural job builder.

Creates 15 z-only constrained sliding-path inputs and a provenance manifest.
DOES NOT perform DFT, contact HPC hosts, create POTCAR, or claim relaxation.
Runs from a fresh checkout: Python stdlib only, no original 48MB db needed.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/"material-screening/source-structures/WSe2_BiDB_row2685_3R_polar_zscan.vasp"
EXPECTED_SHA="34bef6196ab86de3ffff225a93f3caa826ba7e42915817fa061478555cb062f7"
SOURCE_UID="WSe2-1cfbe6183886-2-1_0_0_1--0.67_-0.33"
ARMS={"A":(1/3,-1/3),"B":(-2/3,-1/3),"C":(1/3,2/3)}
TOP=(1,4,5)
STEPS=(0.,.25,.5,.75,1.)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def source_atoms(path=SOURCE):
    path=Path(path)
    if sha(path)!=EXPECTED_SHA:
        raise RuntimeError("Source WSe2 BiDB row2685 POSCAR SHA256 mismatch: abort.")
    s=path.read_text().splitlines()
    if SOURCE_UID not in s[0] or "NOT full relaxation" not in s[0]:
        raise RuntimeError("BiDB row UID or z-scan provenance missing")
    lattice=[list(map(float,s[i].split())) for i in range(2,5)]
    if s[5].split()!=["W","Se"] or s[6].split()!=["2","4"]:
        raise RuntimeError("Species/order not W2Se4")
    atoms=[list(map(float,s[i].split()[:3])) for i in range(8,14)]
    if len(atoms)!=6:raise RuntimeError("Not six atoms")
    if abs(atoms[1][0]-1/3)>1e-7 or abs(atoms[1][1]-2/3)>1e-7:
        raise RuntimeError("Missing expected polar W-on-Se stacking")
    return lattice,atoms

def write_poscar(path,cell,atoms,title,z_only=True):
    header=[title,"1.0",*[" ".join(f"{x:.12f}" for x in v) for v in cell],
            "W Se","2 4"]
    header.extend(["Selective dynamics","Direct"] if z_only else ["Direct"])
    flags=" F F T" if z_only else ""
    lines=header+[" ".join(f"{x%1:.12f}" for x in p)+flags for p in atoms]
    path.write_text("\n".join(lines)+"\n")

def incar(center_z):
    return f"""SYSTEM = WSe2 BiDB2685 pilot z-only structural relaxation; UNRUN
PREC = Accurate
ENCUT = 520
GGA = PE
IVDW = 12
LASPH = .TRUE.
LREAL = .FALSE.
ADDGRID = .TRUE.
ISPIN = 1
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1E-7
EDIFFG = -0.01
IBRION = 2
NSW = 160
ISIF = 2
ISYM = 0
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL = 0.5 0.5 {center_z:.9f}
LWAVE = .FALSE.
LCHARG = .FALSE.
# Local PBE-PAW W Se POTCAR required. ENCUT >= 1.3*max(ENMAX);
# convergence, PBE-D3 damping, and vacuum/dipole must be checked.
"""

def build(out,scope="all"):
    if scope not in ("all","pilot"):raise ValueError("scope must be all/pilot")
    out=Path(out)
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Refusing to overwrite nonempty target directory")
    out.mkdir(parents=True,exist_ok=True)
    cell,orig=source_atoms()
    center=(min(p[2] for p in orig)+max(p[2] for p in orig))/2
    a=math.hypot(*cell[0][:2])
    selection={("A",0.),("A",.25),("A",.5)}
    work=[]
    for arm,delta in ARMS.items():
        for s in STEPS:
            name=f"path_{arm}_s{s:.2f}"
            atoms=[p.copy() for p in orig]
            for i in TOP:
                atoms[i][0]=(atoms[i][0]+s*delta[0])%1
                atoms[i][1]=(atoms[i][1]+s*delta[1])%1
            job=out/name
            job.mkdir()
            write_poscar(job/"POSCAR",cell,atoms,
                f"WSe2 row2685 AB to derived BA: {arm} s={s:.2f}; z only")
            (job/"INCAR").write_text(incar(center))
            (job/"KPOINTS").write_text(
                "15x15x1 Gamma initial mesh, not production converged\n0\nGamma\n15 15 1\n0 0 0\n")
            work.append({"name":name,"arm":arm,"s":s,
                "selected_for_pilot":(arm,s) in selection,
                "status":"NOT_RUN","kind":"constrained_z_only",
                "POSCAR_sha256":sha(job/"POSCAR"),
                "INCAR_sha256":sha(job/"INCAR"),
                "KPOINTS_sha256":sha(job/"KPOINTS")})
    manifest={
        "stage":"8g","scientific_status":"GENERATED_INPUTS_ONLY",
        "source":"original BiDB ASE row2685 z-scan, NOT fully ionically relaxed",
        "source_bilayer_uid":SOURCE_UID,"source_POSCAR_sha256":EXPECTED_SHA,
        "source_bidb_db_sha256":"afe1fe7711189be5b0e189b0c101d6488a7578f4f9ba4e1bb0164534f4255984",
        "a_A":a,"shortest_slide_A":a/math.sqrt(3),
        "arms_fractional":ARMS,"max_initial_pilot_jobs":3,
        "jobs":work,"energies_eV":None,
        "software_status":"not_checked","hpc_status":"NOT_SUBMITTED",
        "caveat":"z-only constrained geometry, no NEB barrier, no SOC gap, no NMEE",
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True)
    p.add_argument("--scope",choices=("all","pilot"),default="all",
                   help="all generates all 15; pilot still generates 15 with 3 flagged")
    args=p.parse_args()
    x=build(args.out,args.scope)
    print(json.dumps({"jobs_prepared":len(x["jobs"]),
        "pilot":[j["name"] for j in x["jobs"] if j["selected_for_pilot"]],
        "status":"NOT_RUN"},indent=2))
