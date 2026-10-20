#!/usr/bin/env python3
"""Stage8g: local-on-HPC guarded PBS pilot. No SSH and no credentials.

Run on the Japan HPC after local harness copies generated inputs. --submit
is the ONLY route to qsub, at most three jobs; reruns skip recorded IDs.
No POTCAR is ever committed. Outputs are NOT NEB/NMEE evidence.
"""
import argparse, hashlib, json, re, shutil, subprocess
from pathlib import Path

PILOT=("path_A_s0.00","path_A_s0.25","path_A_s0.50")

def checksum(path):
    h=hashlib.sha256()
    with open(path,"rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def load_inputs(where):
    root=Path(where).resolve()
    data=json.loads((root/"manifest.json").read_text())
    chosen=[x for x in data["jobs"] if x["selected_for_pilot"]]
    if data["stage"]!="8g" or tuple(x["name"] for x in chosen)!=PILOT:
        raise RuntimeError("Pilot manifest altered")
    for job in chosen:
        folder=root/job["name"]
        for name in ("POSCAR","INCAR","KPOINTS"):
            if checksum(folder/name)!=job[name+"_sha256"]:
                raise RuntimeError("Scientific input changed: "+str(folder/name))
        if "F F T" not in (folder/"POSCAR").read_text():
            raise RuntimeError("Pilot in-plane degrees not fixed")
    return root,chosen

def potcar_gate(file):
    s=Path(file).read_bytes().decode("latin1","replace")
    titles=re.findall(r"(?m)^\s*TITEL\s*=\s*(.+)$",s)
    enmax=[float(x) for x in re.findall(r"\bENMAX\s*=\s*([\d.]+)\s*;",s)]
    if len(titles)!=2 or len(enmax)<2:
        raise RuntimeError("Expected exactly two local PBE-PAW blocks")
    names=[]
    for title in titles:
        if "PAW_PBE" not in title:raise RuntimeError("Non-PBE PAW potential")
        match=re.search(r"\b(W|Se)(?:_sv|_pv)?\b",title)
        names.append(match.group(1) if match else "OTHER")
    if names!=["W","Se"]:raise RuntimeError("POTCAR must be W then Se")
    if 520<1.3*max(enmax)-1e-6:
        raise RuntimeError("ENCUT=520 insufficient for local POTCAR ENMAX")
    return {"sha256":checksum(file),"species":names,"max_ENMAX_eV":max(enmax)}

def load_site(file):
    cfg=json.loads(Path(file).read_text())
    required=("scheduler","queue","ncpus","mem_gb","walltime",
              "vasp_binary","mpi_command","license_confirmed","validated")
    if any(k not in cfg for k in required):raise RuntimeError("Incomplete site profile")
    if cfg["scheduler"]!="pbs" or cfg["validated"] is not True:
        raise RuntimeError("Site not independently verified as PBS")
    if not cfg["license_confirmed"]:
        raise RuntimeError("VASP authorization not verified")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",cfg["queue"]):
        raise RuntimeError("Invalid queue name")
    if not 1<=int(cfg["ncpus"])<=64 or not 2<=int(cfg["mem_gb"])<=256:
        raise RuntimeError("Pilot resources outside guard bounds")
    if not re.fullmatch(r"\d{1,3}:\d\d:\d\d",cfg["walltime"]):
        raise RuntimeError("Invalid walltime")
    binary=Path(cfg["vasp_binary"])
    if not binary.is_absolute() or not binary.is_file():
        raise RuntimeError("VASP binary must exist at verified HPC path")
    if cfg["mpi_command"] not in ("mpirun","mpiexec"):
        raise RuntimeError("Only validated MPI launcher allowed")
    if shutil.which("qsub") is None or shutil.which("qstat") is None:
        raise RuntimeError("Cannot find qsub/qstat on this host")
    bootstrap=cfg.get("bootstrap_script")
    if bootstrap and (not Path(bootstrap).is_absolute() or not Path(bootstrap).is_file()):
        raise RuntimeError("Invalid optional HPC environment bootstrap")
    return cfg

def pbs_script(cfg,name):
    boot=("source "+cfg["bootstrap_script"]+"\n") if cfg.get("bootstrap_script") else ""
    return (f"#!/bin/bash\n#PBS -N {name.replace('.','_')[:15]}\n"
      f"#PBS -q {cfg['queue']}\n"
      f"#PBS -l select=1:ncpus={cfg['ncpus']}:mem={cfg['mem_gb']}gb\n"
      f"#PBS -l walltime={cfg['walltime']}\n#PBS -j oe\n"
      "set -euo pipefail\ncd \"$PBS_O_WORKDIR\"\n"+boot+
      f"{cfg['mpi_command']} -np {cfg['ncpus']} {cfg['vasp_binary']} "
      "> vasp.stdout 2> vasp.stderr\n")

def dispatch(root,site,execute=False):
    root,selected=load_inputs(root)
    cfg=load_site(site)
    potcars=[]
    for j in selected:
        local=root/j["name"]/"POTCAR"
        if not local.exists():raise RuntimeError("Missing LICENSED LOCAL POTCAR: "+str(local))
        potcars.append(potcar_gate(local))
    if len({p["sha256"] for p in potcars})!=1:
        raise RuntimeError("Different pseudopotentials across pilot jobs")
    statefile=root/"pbs_ids.json"
    ids=json.loads(statefile.read_text()) if statefile.exists() else {}
    for job in PILOT:
        file=root/job/"submit.pbs"
        content=pbs_script(cfg,job)
        if file.exists() and file.read_text()!=content:
            raise RuntimeError("An existing PBS script would be overwritten")
        file.write_text(content)
    if not execute:
        return {"status":"DRY_RUN_VALIDATED","pilot":PILOT,
                "queue":cfg["queue"],"ncpus":cfg["ncpus"],
                "potcar_sha256":potcars[0]["sha256"],"existing_job_ids":ids}
    for job in PILOT:
        if job in ids:continue
        result=subprocess.run(["qsub","submit.pbs"],cwd=root/job,
                              text=True,capture_output=True,timeout=30)
        if result.returncode:
            raise RuntimeError("qsub failed "+job+": "+result.stderr[-300:])
        job_id=result.stdout.strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+",job_id):
            raise RuntimeError("Unexpected qsub response; check scheduler for duplicates")
        ids[job]=job_id
        statefile.write_text(json.dumps(ids,indent=2)+"\n")
    return {"status":"SUBMITTED_IDEMPOTENT","pbs_ids":ids,
            "note":"Submission does not imply convergence or completed DFT"}

def read_poscar(file):
    lines=Path(file).read_text().splitlines()
    scale=float(lines[1].strip())
    lat=[[float(v)*scale for v in lines[i].split()[:3]] for i in range(2,5)]
    counts=[int(x) for x in lines[6].split()]
    n=sum(counts)
    coord_line=7
    if lines[coord_line].lower().startswith("s"):coord_line+=1
    if not lines[coord_line].strip().lower().startswith("d"):
        raise RuntimeError("Expected Direct POSCAR/CONTCAR coordinates")
    frac=[[float(v) for v in lines[i].split()[:3]]
          for i in range(coord_line+1,coord_line+1+n)]
    if n!=6 or lines[5].split()!=["W","Se"] or counts!=[2,4]:
        raise RuntimeError("Final structure is not expected six-atom W2Se4")
    return lat,frac

def review_output(folder):
    out=folder/"OUTCAR"
    if not out.exists():return {"status":"NOT_RUN_OR_OUTPUT_NOT_COLLECTED"}
    content=out.read_text(errors="replace")
    vals=re.findall(r"free\s+energy\s+TOTEN\s*=\s*([-+\d.Ee]+)\s+eV",content)
    if not vals:return {"status":"INVALID_NO_TOTEN"}
    if not (folder/"CONTCAR").exists():return {"status":"INCOMPLETE_NO_CONTCAR"}
    rec={"final_TOTEN_eV":float(vals[-1]),
         "OUTCAR_sha256":checksum(out),"CONTCAR_sha256":checksum(folder/"CONTCAR")}
    try:
        old_lat,old=read_poscar(folder/"POSCAR")
        new_lat,new=read_poscar(folder/"CONTCAR")
        cell_error=max(abs(a-b) for x,y in zip(old_lat,new_lat)
                       for a,b in zip(x,y))
        dx=max(abs(((v[0]-u[0]+0.5)%1-0.5)*old_lat[0][0]+
                   ((v[1]-u[1]+0.5)%1-0.5)*old_lat[1][0])
               for u,v in zip(old,new))
        dy=max(abs(((v[1]-u[1]+0.5)%1-0.5)*old_lat[1][1])
               for u,v in zip(old,new))
        rec.update({"fixed_cell_error_A":cell_error,
                    "max_lateral_xy_drift_A":max(dx,dy)})
        headings=list(re.finditer(r"POSITION\s+TOTAL-FORCE\s+\(eV/Angst\)",content))
        if not headings:raise RuntimeError("Missing final force table")
        forces=[]
        for line in content[headings[-1].end():].splitlines():
            fields=line.split()
            if len(fields)==6:
                try:forces.append([float(x) for x in fields])
                except ValueError:continue
            if len(forces)==6:break
        if len(forces)!=6:raise RuntimeError("Final ionic force table incomplete")
        fz=max(abs(f[5]) for f in forces)
        rec["max_abs_Fz_eV_A"]=fz
        rec["vasp_ionic_accuracy_flag"]=("reached required accuracy" in content)
        rec["status"]=("CONVERGED_REQUIRES_SCIENTIFIC_REVIEW"
            if fz<=.012 and rec["vasp_ionic_accuracy_flag"]
            and max(dx,dy)<=.002 and cell_error<=.0001
            else "NOT_CONVERGED_OR_FAILED_CONSTRAINTS")
    except (ValueError,IndexError,RuntimeError) as exc:
        rec["status"]="INVALID_GEOMETRY_OR_FORCES"
        rec["reason"]=str(exc)
    return rec

def audit(root):
    root,selected=load_inputs(root)
    ids_file=root/"pbs_ids.json"
    ids=json.loads(ids_file.read_text()) if ids_file.exists() else {}
    rows=[]
    for j in selected:
        rows.append({"job":j["name"],"pbs_id":ids.get(j["name"]),
                     **review_output(root/j["name"])})
    return {"scope":"initial z-only relaxation, NOT NEB barrier or NMEE",
            "rows":rows,"no_barrier_claim":True,"no_NMEE_claim":True}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--jobs",required=True)
    p.add_argument("--site",help="Local HPC profile (not committed)")
    p.add_argument("--submit",action="store_true")
    p.add_argument("--audit",action="store_true")
    args=p.parse_args()
    if args.audit:ans=audit(args.jobs)
    else:
        if not args.site:p.error("--site required for PBS preflight")
        ans=dispatch(args.jobs,args.site,args.submit)
    print(json.dumps(ans,indent=2))
