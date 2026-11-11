#!/usr/bin/env python3
"""Stage8p: reproducible, zero-new-job real 96-band AB/GLIDE comparison.

Uses existing Japan IMR Stage8m EIG/MMN/SPN/WIN only. It fixes Stage8o's
factor-two mismatch (-mu_B*B*sigma_z vs -2mu_B*B*sigma_z), source 1/2
normalization, and wrong identification of <u|du/dk> with dP/dk.
Uses local points actually ON THE 21x21 k mesh, not mislabeled coordinates.

Output contains only input hashes and sanitized numbers. It does NOT
establish full-BZ Gamma or physical layer-odd response.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from stage8m_plaq import read_eig, read_mmn_spn, link_matrix, DeterminantPlaquette
from stage8m_response import (
    read_win_cell,read_win_kpoints,read_spn_formatted,WannierProvider)
from stage8p_projector_kubo_core import compare_local_three_way,MU_B

OCC=36
SAMPLES=((8,2),(7,9))
DL=1e-4
DB=.25


def hash_sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for buf in iter(lambda:f.read(8*1024*1024),b""):h.update(buf)
    return h.hexdigest()


def run(d):
    d=Path(d)
    paths={n:d/f"wannier90.{n}" for n in ("eig","mmn","spn","win")}
    for n,path in paths.items():
        if not path.is_file():raise FileNotFoundError(str(path))
    E=read_eig(paths["eig"])
    nb,nk,nn,mats=read_mmn_spn(paths["mmn"])
    S=read_spn_formatted(paths["spn"],nb,nk)
    kpts=read_win_kpoints(paths["win"])
    if (nb,nk)!=(96,441) or E.shape!=(441,96) or kpts.shape[0]!=441:
        raise ValueError("Unexpected 96-band x 441-k source geometry")
    kmesh=np.rint(kpts[:,:2]*21).astype(int)
    mask=np.abs(((kpts[:,:2]*21-kmesh+10.5)%21)-10.5)
    if np.max(mask)>1e-4:raise ValueError("k points do not lie on strict 21x21 lattice")
    idx={}
    for i,(a,b) in enumerate(kmesh):
        key=(int(a%21),int(b%21))
        if key in idx:raise ValueError("Duplicate folded k")
        idx[key]=i
    if len(idx)!=441:raise ValueError("Incomplete k mesh")
    neigh={(a,s):{} for a in (0,1) for s in (-1,1)}
    for (u,v),i in idx.items():
        neigh[(0,+1)][i]=idx[((u+1)%21,v)]
        neigh[(0,-1)][i]=idx[((u-1)%21,v)]
        neigh[(1,+1)][i]=idx[(u,(v+1)%21)]
        neigh[(1,-1)][i]=idx[(u,(v-1)%21)]
    cell=read_win_cell(paths["win"])
    a1,a2=cell[0,:2],cell[1,:2]
    area=float(np.linalg.det(np.array([a1,a2])))
    b1=2*np.pi*np.array([a2[1],-a2[0]])/area
    b2=2*np.pi*np.array([-a1[1],a1[0]])/area
    dk=(float(np.linalg.norm(b1))/21,float(np.linalg.norm(b2))/21)
    rows=[]
    for grid in SAMPLES:
        ik=idx[grid]
        for src in (0,1):
            srcmat=S[:,:,:,src]
            provider=WannierProvider(E,srcmat,S[:,:,:,2],OCC,nb)
            for ax in (0,1):
                jp=neigh[(ax,+1)][ik];jm=neigh[(ax,-1)][ik]
                Mp=link_matrix(mats,ik,jp,nb);Mm=link_matrix(mats,ik,jm,nb)
                if Mp.shape!=(nb,nb) or Mm.shape!=(nb,nb):
                    raise ValueError("Missing/mismatched MMN bond")
                tuples=[
                    (E[n],S[n,:,:,src],S[n,:,:,2])
                    for n in (ik,jp,jm)
                ]
                vals={}
                maxleak=0
                maxprojector=0
                for B in (DB,-DB):
                    r=compare_local_three_way(tuples[0],tuples[1],tuples[2],
                        (Mp,Mm),B=B,occ=OCC,dk=dk[ax],dl=DL)
                    vals[B]=r
                    maxleak=max(maxleak,r["MMN_subspace_leakage"]["max_link_nonunitarity_fro"])
                    maxprojector=max(maxprojector,r["MMN_subspace_leakage"]["max_neighbour_projector_idempotency_error"])
                Gres=(vals[DB]["F_resolvent"]-vals[-DB]["F_resolvent"])/(2*DB)
                Gproj=(vals[DB]["F_projector_fd"]-vals[-DB]["F_projector_fd"])/(2*DB)
                eng=DeterminantPlaquette(provider.frames,mode="mmn",links=mats,
                    ik_of=lambda k,s,a=ax:neigh[(a,1 if s>0 else -1)][k],
                    dk=dk[ax],dl=DL,axis=ax,full_nb=nb)
                Gw=eng.Gamma_k(ik,0,0,DB)
                pplus=eng.plaquette_phase(ik,0,DB)
                pminus=eng.plaquette_phase(ik,0,-DB)
                rows.append({
                    "sample_grid_indices":list(grid),
                    "actual_k_frac":list(map(float,kpts[ik][:2])),
                    "source":"sigma_x/2" if src==0 else "sigma_y/2",
                    "axis_oblique_reciprocal":ax,
                    "k_index_0based":ik,
                    "occ_bands":OCC,
                    "gap_center_eV":vals[DB]["occ_unocc_gap_eV"],
                    "gamma_resolvent":Gres,
                    "gamma_projector":Gproj,
                    "gamma_original_wilson":Gw,
                    "relative_resolvent_projector_difference":abs(Gres-Gproj)/max(abs(Gres),abs(Gproj),1e-30),
                    "relative_resolvent_wilson_difference":abs(Gres-Gw)/max(abs(Gres),abs(Gw),1e-30),
                    "max_link_nonunitarity_fro":maxleak,
                    "max_neighbor_projector_idempotency_error":maxprojector,
                    "raw_signed_loop_angle_difference_rad":(pplus-pminus)*dk[ax]*DL,
                    "raw_phase_diff_to_machine_eps":abs((pplus-pminus)*dk[ax]*DL)/np.finfo(float).eps
                })
    return {
        "status":"ZERO_JOB_REAL_LOCAL_KUBO_PROJECTOR_WILSON_COMPARISON_NOT_FULL_BZ_GAMMA",
        "new_PBS_jobs":0,
        "input_sha256":{n:hash_sha(p) for n,p in paths.items()},
        "sources":"lambda*sigma_{x/y}/2; Zeeman -mu_B*B*sigma_z, no second g factor",
        "reciprocal_derivative":"oblique b1/b2 unit direction; NOT cartesian tensor without full transform",
        "k_points":len(kpts),"bands":nb,
        "limitations":["Finite 96-band truncated MMN links may be nonunitary",
                       "All three estimates derive from same finite-band VASP files",
                       "No layer-odd spin projector or true MMN symmetry sewing",
                       "No full-BZ integration or accepted WSe2 Gamma"],
        "rows":rows
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--dir",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    out=Path(args.out)
    if out.exists():raise FileExistsError("Refusing to overwrite prior results")
    ans=run(args.dir)
    out.write_text(json.dumps(ans,indent=2)+"\n")
    print(json.dumps({"status":ans["status"],
          "number_of_local_cross_checks":len(ans["rows"]),
          "max_relative_cross_method_discrepancy":max(x["relative_resolvent_wilson_difference"] for x in ans["rows"])},
          indent=2))
