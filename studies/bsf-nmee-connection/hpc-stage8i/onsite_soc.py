#!/usr/bin/env python3
"""Stage8i: onsite guarded two-job SOC PBS driver, plus read-only status/audit.

No network, no SSH credentials, no licensed potentials stored in repo.
Local DeepSeek harness is responsible for site discovery, staging and monitor.
"""
import argparse,hashlib,json,re,shutil,subprocess
from pathlib import Path

NAMES=("soc_AB_s0.00","soc_GLIDE_s0.50")
EXPECTED_POTCAR_SHA="7ef8b236d47cb52384021b0d5255b9dd2cd6fb636f1b7e0a093eff353b3a2a5c"

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):h.update(data)
    return h.hexdigest()

def verify_inputs(root):
    root=Path(root).resolve()
    data=json.loads((root/"manifest.json").read_text())
    if data.get("stage")!="8i" or tuple(data.get("jobs",{}))!=NAMES:
        raise RuntimeError("Stage8i must contain exactly two jobs")
    for name in NAMES:
        folder=root/name
        d=data["jobs"][name]
        if d["POSCAR_sha256"]!=d["stage8h_original_contcar_sha256"]:
            raise RuntimeError("Stage8h source structure changed for "+name)
        for key in ("POSCAR","INCAR","KPOINTS"):
            if digest(folder/key)!=d[key+"_sha256"]:
                raise RuntimeError("Hash mismatch "+name+"/"+key)
        incar=(folder/"INCAR").read_text()
        for tag in ("LSORBIT = .TRUE.","LNONCOLLINEAR = .TRUE.",
                    "NSW = 0","ISYM = 0","NBANDS = 96"):
            if tag not in incar:raise RuntimeError("Missing spinor SCF tag "+tag)
        if re.search(r"(?m)^\s*ISPIN\s*=",incar):
            raise RuntimeError("ISPIN in noncollinear spinor INCAR is forbidden")
    return root,data

def potcar_gate(path):
    content=Path(path).read_bytes().decode("latin1","replace")
    title=re.findall(r"(?m)^\s*TITEL\s*=\s*(.+)$",content)
    enmax=[float(v) for v in re.findall(r"\bENMAX\s*=\s*([\d.]+)\s*;",content)]
    if len(title)!=2 or len(enmax)<2:raise RuntimeError("Missing PBE PAW W/Se blocks")
    specie=[]
    for t in title:
        if "PAW_PBE" not in t:raise RuntimeError("Wrong PAW/XC family")
        m=re.search(r"\b(W|Se)(?:_sv|_pv)?\b",t)
        specie.append(m.group(1) if m else "invalid")
    if specie!=["W","Se"]:raise RuntimeError("POTCAR species must be W then Se")
    if max(enmax)*1.3>520+1e-8:raise RuntimeError("ENCUT 520 insufficient")
    d=digest(path)
    if d!=EXPECTED_POTCAR_SHA:
        raise RuntimeError("POTCAR differs from verified Stage8g W/Se reference; resolve before running")
    return {"sha256":d,"species":specie,"max_enmax_eV":max(enmax)}

def site_gate(file):
    cfg=json.loads(Path(file).read_text())
    required=("scheduler","validated","license_confirmed","soc_binary_verified",
              "queue","ncpus","mem_gb","walltime","vasp_binary","mpi_command")
    if any(k not in cfg for k in required):raise RuntimeError("Missing site keys")
    if cfg["scheduler"]!="pbs" or not cfg["validated"] or not cfg["license_confirmed"]:
        raise RuntimeError("Unvalidated PBS license/site")
    if cfg["soc_binary_verified"] is not True:
        raise RuntimeError("VASP executable not verified for noncollinear LSORBIT")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",cfg["queue"]):
        raise RuntimeError("Unsafe queue")
    if not 1<=int(cfg["ncpus"])<=64 or not 2<=int(cfg["mem_gb"])<=256:
        raise RuntimeError("Unsafe resource demand")
    if not re.fullmatch(r"\d{1,3}:\d\d:\d\d",cfg["walltime"]):
        raise RuntimeError("Bad walltime")
    binary=Path(cfg["vasp_binary"])
    if not binary.is_absolute() or not binary.is_file():
        raise RuntimeError("Actual onsite SOC binary required")
    if "ncl" not in binary.name and not cfg.get("explicitly_tested_soc_on_std_binary",False):
        raise RuntimeError("Need VASP noncollinear/SOC capable binary (normally vasp_ncl)")
    if cfg["mpi_command"] not in ("mpirun","mpiexec"):
        raise RuntimeError("MPI command not verified")
    if shutil.which("qsub") is None or shutil.which("qstat") is None:
        raise RuntimeError("PBS not available onsite")
    boot=cfg.get("bootstrap_script")
    if boot and (not Path(boot).is_absolute() or not Path(boot).is_file()):
        raise RuntimeError("Missing actual bootstrap script")
    return cfg

def script(cfg,name):
    boot=("source "+cfg["bootstrap_script"]+"\n") if cfg.get("bootstrap_script") else ""
    return (f"#!/bin/bash\n#PBS -N {name[:15]}\n#PBS -q {cfg['queue']}\n"
      f"#PBS -l select=1:ncpus={cfg['ncpus']}:mem={cfg['mem_gb']}gb\n"
      f"#PBS -l walltime={cfg['walltime']}\n#PBS -j oe\n"
      "set -euo pipefail\ncd \"$PBS_O_WORKDIR\"\n"+boot+
      f"{cfg['mpi_command']} -np {int(cfg['ncpus'])} {cfg['vasp_binary']} "
      "> vasp.stdout 2> vasp.stderr\n")

def launch(root,site,submit=False):
    root,manifest=verify_inputs(root)
    cfg=site_gate(site)
    infos=[]
    for name in NAMES:
        pp=root/name/"POTCAR"
        if not pp.exists():raise RuntimeError("Missing authorized onsite POTCAR: "+name)
        infos.append(potcar_gate(pp))
    if len({x["sha256"] for x in infos})!=1:
        raise RuntimeError("Different PBE PAW files across SOC cases")
    state_file=root/"soc_pbs_ids.json"
    state=json.loads(state_file.read_text()) if state_file.exists() else {}
    for name in NAMES:
        job=root/name/"submit_soc.pbs"
        generated=script(cfg,name)
        if job.exists() and job.read_text()!=generated:
            raise RuntimeError("Existing PBS script differs; stop for explicit review")
        job.write_text(generated)
    if not submit:
        return {"status":"VALIDATED_DRY_RUN","jobs":list(NAMES),
                "queue":cfg["queue"],"vasp_binary":cfg["vasp_binary"],
                "potcar_sha256":infos[0]["sha256"],
                "existing_pbs_ids":state}
    for name in NAMES:
        if name in state:continue
        # Fail closed on interrupted qsub: agent must resolve scheduler state
        # before retrying a potentially accepted submission.
        state[name]="SUBMITTING_CONFIRM_SCHEDULER_IF_INTERRUPTED"
        state_file.write_text(json.dumps(state,indent=2)+"\n")
        result=subprocess.run(["qsub","submit_soc.pbs"],cwd=root/name,
            text=True,capture_output=True,timeout=40)
        if result.returncode:
            raise RuntimeError("qsub failed; inspect scheduler before clearing SUBMITTING state: "+name)
        jobid=result.stdout.strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+",jobid):
            raise RuntimeError("Unexpected PBS response; inspect queue before retry")
        state[name]=jobid
        state_file.write_text(json.dumps(state,indent=2)+"\n")
    return {"status":"SUBMITTED_OR_RESUMED","pbs_ids":state,
            "note":"A PBS ID is neither completed SCF nor finite mixed Berry response"}

def status(root):
    root,_=verify_inputs(root)
    s=root/"soc_pbs_ids.json"
    state=json.loads(s.read_text()) if s.is_file() else {}
    rows=[]
    for name in NAMES:
        jobid=state.get(name)
        result={"job":name,"pbs_id":jobid}
        if not jobid:
            result["state"]="NOT_SUBMITTED"
        elif jobid.startswith("SUBMITTING_"):
            result["state"]="UNKNOWN_SUBMISSION_CHECK_SCHEDULER"
        else:
            p=subprocess.run(["qstat","-f",jobid],text=True,
                             capture_output=True,timeout=25)
            match=re.search(r"\bjob_state\s*=\s*([A-Z])",p.stdout)
            result["state"]=match.group(1) if p.returncode==0 and match else "UNKNOWN_OR_PURGED"
        rows.append(result)
    return {"stage":"8i","rows":rows,"no_soc_gap_or_NMEE_inferred_from_scheduler":True}

def inspect(root):
    root,manifest=verify_inputs(root)
    result=[]
    for name in NAMES:
        folder=root/name
        out=folder/"OUTCAR"
        row={"job":name,"status":"NO_OUTPUT_OR_NOT_FINISHED","total_static_NMEE":None}
        if out.exists():
            s=out.read_text(errors="replace")
            en=re.findall(r"free\s+energy\s+TOTEN\s*=\s*([-+\d.Ee]+)\s+eV",s)
            if en:
                row["TOTEN_eV"]=float(en[-1])
                row["status"]="SCF_OUTPUT_FOUND_FURTHER_CHECKS_REQUIRED"
                row["OUTCAR_sha256"]=digest(out)
            row["EIGENVAL_present"]=(folder/"EIGENVAL").exists()
            row["PROCAR_present"]=(folder/"PROCAR").exists()
            row["WAVECAR_present"]=(folder/"WAVECAR").exists()
            if (folder/"EIGENVAL").exists():
                from eigenval_gate import analyze_eigenval
                row["band_analysis"]=analyze_eigenval(folder/"EIGENVAL",21*21)
        result.append(row)
    return {"scope":"two 21x21x1 SOC SCF candidate electronic feasibility gates",
            "jobs":result,
            "spin_mixing_matrix_element_status":"NOT_PROVEN_BY_PROCAR_DIAGONAL_EXPECTATIONS",
            "material_static_NMEE":None,
            "claim_of_continuous_BZ_gap":False}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--jobs",required=True);p.add_argument("--site")
    p.add_argument("--submit",action="store_true")
    p.add_argument("--status",action="store_true")
    p.add_argument("--audit",action="store_true")
    a=p.parse_args()
    if a.status: result=status(a.jobs)
    elif a.audit:result=inspect(a.jobs)
    else:
        if not a.site:p.error("Local site JSON needed for dry run/submit")
        result=launch(a.jobs,a.site,submit=a.submit)
    print(json.dumps(result,indent=2))
