#!/usr/bin/env python3
"""Stage 8i: immutable two-job SOC preflight built from ACTUAL Stage8h CONTCARs.

Run where authenticated Stage8g files exist, not from fabricated structures.
No SSH/POTCAR/qsub. Safe standard-library-only builder; all energies NOT_RUN.
"""
import argparse, hashlib, json, re, shutil
from pathlib import Path

EXPECTED={
    "soc_AB_s0.00": {"stage8g":"path_A_s0.00",
                     "sha256":"60ee4c590324dea8b7e3fa99780dc496717f12bc29f3a4f7deffe423d3e57992",
                     "operations":6,"total_tensor_dim":1},
    "soc_GLIDE_s0.50":{"stage8g":"path_A_s0.50",
                      "sha256":"fe5eb6f6071f6a8f9a1dd76d4603b1dad694f1fbaf5d854bb49107ce79b9fb5b",
                      "operations":4,"total_tensor_dim":0},
}

def digest(path):
    sha=hashlib.sha256()
    with Path(path).open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):sha.update(piece)
    return sha.hexdigest()

def parse_poscar(path):
    lines=Path(path).read_text().splitlines()
    if len(lines)<14:raise ValueError("CONTCAR truncated")
    scale=float(lines[1].split()[0])
    if scale<=0:raise ValueError("Only positive POSCAR scale supported")
    cell=[[float(w)*scale for w in lines[j].split()[:3]] for j in range(2,5)]
    if lines[5].split()!=["W","Se"] or lines[6].split()!=["2","4"]:
        raise ValueError("Expected source W2Se4 with W then Se")
    start=7
    if lines[start].lstrip().lower().startswith("s"):start+=1
    tag=lines[start].strip().lower()
    if not tag.startswith("d"):raise ValueError("Expected Direct CONTCAR")
    sites=[[float(x) for x in row.split()[:3]]
           for row in lines[start+1:start+7]]
    if len(sites)!=6 or any(len(site)!=3 for site in sites):
        raise ValueError("Not six direct-coordinate sites")
    c=cell[2][2]
    return cell,sites,sum(site[2] for site in sites)/6%1,lines

def incar(dipole_z):
    return f"""SYSTEM = WSe2 Stage8i nonmagnetic spinor-SOC feasibility SCF (NOT NMEE)
PREC = Accurate
GGA = PE
ENCUT = 520
IVDW = 12
LASPH = .TRUE.
LREAL = .FALSE.
ADDGRID = .TRUE.
ALGO = Normal
NELM = 140
EDIFF = 1E-7
ISTART = 0
ICHARG = 2
ISMEAR = 0
SIGMA = 0.05
NSW = 0
IBRION = -1
ISYM = 0
LNONCOLLINEAR = .TRUE.
LSORBIT = .TRUE.
SAXIS = 0 0 1
MAGMOM = 18*0
NBANDS = 96
LORBIT = 11
LWAVE = .TRUE.
LCHARG = .TRUE.
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL = 0.5 0.5 {dipole_z:.9f}
# Do not add ISPIN=2; use an ACTUAL VASP noncollinear-capable executable.
# W/Se licensed local PBE PAW POTCAR; no POTCAR is distributed by this repo.
# NBANDS=96 is an initial transition window, not a Wannier convergence claim.
# k=21x21x1 is a sampled 2D mesh, not a mathematically guaranteed full BZ gap.
"""

KPOINTS="21x21x1 Gamma SOC sampled whole 2D BZ\n0\nGamma\n21 21 1\n0 0 0\n"

def build(existing_stage8g_root,out_dir):
    source_root=Path(existing_stage8g_root)
    out=Path(out_dir)
    if out.exists():
        if any(out.iterdir()):
            raise RuntimeError("Refuse overwriting existing SOC work")
    out.mkdir(parents=True,exist_ok=True)
    refs={}
    for name,meta in EXPECTED.items():
        src=source_root/meta["stage8g"]/"CONTCAR"
        if not src.is_file():raise FileNotFoundError("Missing authentic Stage8h CONTCAR: "+str(src))
        got=digest(src)
        if got!=meta["sha256"]:
            raise RuntimeError(f"FAILED_STAGE8H_CONTCAR_SHA: {name} {got}")
        cell,coords,dip,lines=parse_poscar(src)
        folder=out/name
        folder.mkdir()
        shutil.copyfile(src,folder/"POSCAR")
        (folder/"INCAR").write_text(incar(dip))
        (folder/"KPOINTS").write_text(KPOINTS)
        refs[name]={"previous_job":meta["stage8g"],
                    "stage8h_original_contcar_sha256":meta["sha256"],
                    "source_group_operations":meta["operations"],
                    "source_total_static_tensor_dim":meta["total_tensor_dim"],
                    "job_status":"NOT_SUBMITTED",
                    "POSCAR_sha256":digest(folder/"POSCAR"),
                    "INCAR_sha256":digest(folder/"INCAR"),
                    "KPOINTS_sha256":digest(folder/"KPOINTS"),
                    "dipole_frac_z":dip}
    manifest={"stage":"8i","purpose":"two REAL spinor-SOC SCF feasibility jobs, NOT NMEE or NEB",
       "max_new_jobs":2,"source":"sha-verified original Stage8h relaxed WSe2 CONTCAR files",
       "run_protocol":"fixed ionic geometry PBE+D3(BJ) SOC; initial 21x21x1 gamma mesh",
       "analysis_limits":"spin expectation and sampled EIGENVAL gap do not establish interband spin transition matrix elements or finite Gamma",
       "jobs":refs,"all_jobs_submitted":False,"actual_soc_results":None,
       "author_original_prl_reproduction_modified":False}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--stage8g-root",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    result=build(args.stage8g_root,args.out)
    print(json.dumps({"status":"PREPARED_NOT_SUBMITTED",
       "jobs":list(result["jobs"]),"source_authentication":"stage8h original CONTCAR SHA256 verified"},indent=2))
