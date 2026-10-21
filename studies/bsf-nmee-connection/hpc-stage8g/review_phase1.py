#!/usr/bin/env python3
"""Stage8h: RAW-OUTPUT science review of three actual Stage8g pilot jobs.

Requires authorized local copy of real OUTCAR/CONTCAR, never commits them.
Checks whole six-atom z-relaxed point-group operations and layer glide.
This is independent from the local-harness summary: FAIL CLOSED on missing
raw evidence, and never labels a 3-point scan as an NEB barrier or NMEE.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
SYMM=STUDY/"material-screening"/"ideal_stack_symmetry.py"
SPEC=importlib.util.spec_from_file_location("stage8h_symmetry",SYMM)
motif=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(motif)
sys.path.insert(0,str(HERE))
import onsite_pbs

GROUP_EXPECTED={"path_A_s0.00":(6,1,False),
                "path_A_s0.25":(2,2,False),
                "path_A_s0.50":(4,0,True)}

def structure_ops(poscar,tol=1e-3):
    cell,frac=onsite_pbs.read_poscar(poscar)
    lat=np.asarray(cell,float)
    atoms=np.asarray(frac,float)
    if lat.shape!=(3,3) or abs(np.linalg.det(lat))<1:
        raise ValueError("Invalid 3D lattice")
    # Origin-independent z shift: the two W centers set the exchange mirror.
    site=atoms.copy()
    site[:,2]=(atoms@lat)[:,2]
    site[:,2]-=(site[0,2]+site[1,2])/2
    elems=["W","W","Se","Se","Se","Se"]
    refl=np.diag([-1.,1.])
    ops=[]
    for n in range(6):
        angle=n*np.pi/3
        rot=np.array([[np.cos(angle),-np.sin(angle)],
                      [np.sin(angle),np.cos(angle)]])
        for mirrored in (False,True):
            R=rot@(refl if mirrored else np.eye(2))
            for fz in (1,-1):
                if motif.maps(site,elems,lat[:2,:2],
                              site,elems,R,fz,tol=tol):
                    ops.append((n,mirrored,fz,R))
    if not ops:raise RuntimeError("No identity operation; bad raw structure")
    glide=any(n==0 and not m and fz==-1 for n,m,fz,R in ops)
    rank={}
    for staggered in (False,True):
        images=[]
        for i in range(2):
            for j in range(2):
                seed=np.zeros((2,2));seed[i,j]=1
                proj=sum((1 if staggered else fz)*R@seed@R.T
                         for n,m,fz,R in ops)/len(ops)
                images.append(proj.reshape(4))
        rank["layer_odd" if staggered else "total"]=int(
            np.linalg.matrix_rank(np.stack(images),tol=1e-8))
    return {"operations":len(ops),"horizontal_glide":glide,
            "tensor_dimension":rank,"symmetry_tolerance":tol,
            "structure_sha256":onsite_pbs.checksum(poscar)}

def gate(root,tol=1e-3):
    root,selected=onsite_pbs.load_inputs(root)
    cases=[]
    for job in selected:
        p=root/job["name"]
        output=onsite_pbs.review_output(p)
        if output["status"]!="CONVERGED_REQUIRES_SCIENTIFIC_REVIEW":
            return {"status":"BLOCKED_RAW_OUTPUT_INCOMPLETE_OR_IONIC_CHECK_FAILED",
                    "job":job["name"],"reported":output,
                    "no_NEBar_or_NMEE_claim":True}
        info=structure_ops(p/"CONTCAR",tol=tol)
        expected_ops,dim,glide=GROUP_EXPECTED[job["name"]]
        info["symmetry_expected_from_ideal"]={
            "operations":expected_ops,"total_tensor_dimension":dim,
            "horizontal_glide":glide}
        info["symmetry_pass"]=(
            info["operations"]==expected_ops and
            info["tensor_dimension"]["total"]==dim and
            info["horizontal_glide"]==glide)
        cases.append({"job":job["name"],"slide_fraction":job["s"],
                      "raw_output_checks":output,**info})
    energies=[r["raw_output_checks"]["final_TOTEN_eV"] for r in cases]
    delta=[1000*(x-energies[0]) for x in energies]
    symmetry_ok=all(c["symmetry_pass"] for c in cases)
    return {"status":("RAW_GEOMETRY_SYMMETRY_GATE_PASS__MESH_XC_CONVERGENCE_PENDING"
                      if symmetry_ok else
                      "RELAXED_SYMMETRY_NOT_CONFIRMED_STOP_BEFORE_NMEE"),
            "source":"local actual OUTCAR + CONTCAR and SHA-locked manifest",
            "jobs":cases,"E_eV_per_2fu_cell":energies,
            "delta_E_meV_per_2fu_cell":delta,
            "E_second_finite_difference_meV":
                 1000*(energies[2]-2*energies[1]+energies[0]),
            "known_missing":["k-mesh/ENCUT/dispersion/vacuum convergence",
                            "mirror partner independent energy check",
                            "full ionic endpoints and NEB barrier",
                            "SOC full BZ gap and real static NMEE"],
            "NEB_barrier":None,"material_static_NMEE":None}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--jobs",required=True)
    parser.add_argument("--tol",type=float,default=.001)
    parser.add_argument("--write",help="optional LOCAL output JSON; never automatically git add")
    p=parser.parse_args()
    answer=gate(p.jobs,p.tol)
    print(json.dumps(answer,indent=2))
    if p.write:
        dest=Path(p.write)
        if dest.exists():raise RuntimeError("Will not overwrite prior audit")
        dest.write_text(json.dumps(answer,indent=2)+"\n")
