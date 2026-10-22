#!/usr/bin/env python3
"""Stage8l: zero-DFT pseudo-spin INTERBAND vertex pilot from existing WAVECAR.

This is a falsifiable low-cost screen only: no PAW augmentation, velocity,
layer operator, Berry connection, Zeeman derivative, BZ integral, or NMEE.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from stage8k_independent_reader import (
    read_header, kpoint_meta, band_coeffs, split_spinor, find_kindex
)

TARGETS={
 "GAMMA_M_TWO_THIRDS_NOT_KPRIME":(-1/3,1/3,0),
 "NEAR_M_NOT_EXACT_M":(10/21,0,0),
 "GAMMA_K_INTERIOR_NOT_CBM":(3/21,3/21,0),
 "TRUE_K":(1/3,1/3,0),
}
# VASP noncollinear spinor bands use one occupation per band: occupied ~1.
OCCUPIED_BANDS=(33,34,35,36)
EMPTY_BANDS=tuple(range(37,45))
K_TOL=1e-5
ENERGY_TOL_EV=1e-5
OCCUPIED_MIN=0.90
EMPTY_MAX=0.10

def sha256_file(path, block_size=8*1024*1024):
    digest=hashlib.sha256()
    with open(path,"rb") as stream:
        for block in iter(lambda:stream.read(block_size),b""):
            digest.update(block)
    return digest.hexdigest()

def parse_eigenval(path):
    """Strict standard-library EIGENVAL parser for noncollinear VASP rows."""
    lines=Path(path).read_text(errors="strict").splitlines()
    if len(lines)<7: raise ValueError("EIGENVAL is truncated before k-point data")
    header=lines[5].split()
    if len(header)<3: raise ValueError("EIGENVAL header lacks NELECT/NKPTS/NBANDS")
    nelect=float(header[0]); nk=int(header[1]); nb=int(header[2])
    if nk<=0 or nb<=0: raise ValueError("Invalid EIGENVAL k-point/band counts")
    i=6;kpts=[];bands=[]
    for ik in range(nk):
        while i<len(lines) and not lines[i].strip():i+=1
        if i>=len(lines):raise ValueError(f"EIGENVAL ends before k-point {ik}")
        row=lines[i].split();i+=1
        if len(row)<4:raise ValueError(f"Malformed EIGENVAL k-point row {ik}")
        kpts.append([float(x) for x in row[:3]])
        br=[]
        for ib in range(nb):
            while i<len(lines) and not lines[i].strip():i+=1
            if i>=len(lines):raise ValueError(f"EIGENVAL ends at k={ik}, band={ib+1}")
            vals=lines[i].split();i+=1
            if len(vals)<3:raise ValueError(f"Malformed EIGENVAL band row k={ik}, band={ib+1}")
            if int(vals[0])!=ib+1:raise ValueError(f"Unexpected band index at k={ik}: {vals[0]}")
            br.append([float(vals[1]),float(vals[2])])
        bands.append(br)
    return kpts,bands,{"nelect":nelect,"nk":nk,"nb":nb}

def source_vertices(pseudo_valence,pseudo_conduction):
    """Return normalized pseudo <v|sigma/2|c>, using the plane-wave Euclidean metric."""
    vu,vd=pseudo_valence;cu,cd=pseudo_conduction
    if vu.shape!=cu.shape or vd.shape!=cd.shape or vu.shape!=vd.shape:
        raise ValueError("Spinor coefficient dimensions differ; common G basis not established")
    nv=float(np.vdot(vu,vu).real+np.vdot(vd,vd).real)
    nc=float(np.vdot(cu,cu).real+np.vdot(cd,cd).real)
    if min(nv,nc)<=1e-10:raise ValueError("Invalid pseudo spinor norm")
    den=(nv*nc)**0.5
    ud=np.vdot(vu,cd);du=np.vdot(vd,cu)
    sx=(ud+du)/(2*den)
    sy=(-1j*ud+1j*du)/(2*den)
    sz=(np.vdot(vu,cu)-np.vdot(vd,cd))/(2*den)
    return np.array([sx,sy,sz])

def screen(wavecar,eigenval,occupied=OCCUPIED_BANDS,empty=EMPTY_BANDS,
           expected_wavecar_sha256=None,expected_eigenval_sha256=None):
    whash=sha256_file(wavecar);ehash=sha256_file(eigenval)
    if expected_wavecar_sha256 and whash.lower()!=expected_wavecar_sha256.lower():
        raise ValueError("WAVECAR SHA256 mismatch against frozen evidence manifest")
    if expected_eigenval_sha256 and ehash.lower()!=expected_eigenval_sha256.lower():
        raise ValueError("EIGENVAL SHA256 mismatch against frozen evidence manifest")
    h=read_header(wavecar)
    kpts,bands,meta=parse_eigenval(eigenval)
    if h["nk"]!=len(kpts) or h["nb"]!=len(bands[0]):
        raise ValueError("WAVECAR/EIGENVAL k/band counts not equal")
    if max(occupied)>=min(empty) or min(occupied)<1 or max(empty)>h["nb"]:
        raise ValueError("Occupied-empty window invalid")
    karr=np.asarray(kpts,float)
    energies_all=np.asarray([[x[0] for x in row] for row in bands],float)
    occ_all=np.asarray([[x[1] for x in row] for row in bands],float)
    results={}
    for name,coord in TARGETS.items():
        idx,distance=find_kindex(h,coord,karr)
        if distance>K_TOL:
            raise ValueError(f"Requested k not represented exactly in EIGENVAL grid: {name}, nearest distance {distance}")
        nplane,wk,enocc=kpoint_meta(wavecar,h,idx)
        wk=np.asarray(wk,float)
        kdelta=np.max(np.abs(((wk-karr[idx]+0.5)%1)-0.5))
        if kdelta>K_TOL:
            raise ValueError(f"WAVECAR/EIGENVAL k-coordinate mismatch at {name}: {kdelta}")
        eig=energies_all[idx];occ=occ_all[idx]
        mismatch=float(np.max(np.abs(enocc[:,0]-eig)))
        if mismatch>ENERGY_TOL_EV:
            raise ValueError(f"WAVECAR/EIGENVAL energy mismatch at {name}: {mismatch}")
        bad_occ=[b for b in occupied if occ[b-1]<OCCUPIED_MIN]
        bad_empty=[b for b in empty if occ[b-1]>EMPTY_MAX]
        if bad_occ or bad_empty:
            raise ValueError(f"Selected window is not fully insulating at {name}: "
                             f"not-occupied={bad_occ}, not-empty={bad_empty}")
        components={}
        for band in sorted(set(occupied).union(empty)):
            coeffs,_=band_coeffs(wavecar,h,idx,band-1)
            up,down,_=split_spinor(coeffs,nplane,"grouped")
            components[band]=(up,down)
        pairs=[]
        for v in occupied:
            for c in empty:
                vec=source_vertices(components[v],components[c])
                pairs.append({
                    "occupied_band_1based":v,"empty_band_1based":c,
                    "delta_E_eV":round(float(eig[c-1]-eig[v-1]),7),
                    "pseudo_abs_Sx":round(float(abs(vec[0])),8),
                    "pseudo_abs_Sy":round(float(abs(vec[1])),8),
                    "pseudo_abs_Sz":round(float(abs(vec[2])),8),
                    "complex_Sx":[round(float(vec[0].real),8),round(float(vec[0].imag),8)]
                })
        # Full-block Frobenius norm is invariant to unitary rotations within
        # the selected occupied and empty subspaces; individual pair entries
        # are not invariant within degenerate subspaces.
        magnitude=float(np.sqrt(sum(p["pseudo_abs_Sx"]**2+p["pseudo_abs_Sy"]**2
                                    for p in pairs)))
        results[name]={
            "requested_fractional_k":list(coord),
            "sampled_fractional_k":[float(a) for a in kpts[idx]],
            "wavecar_fractional_k":[float(a) for a in wk],
            "k_index_0based":idx,"k_distance_max_mod1":float(distance),
            "wavecar_eigenval_k_distance_max_mod1":float(kdelta),
            "source_energy_roundtrip_max_diff_eV":mismatch,
            "selected_occupied_occupations_min":float(min(occ[b-1] for b in occupied)),
            "selected_empty_occupations_max":float(max(occ[b-1] for b in empty)),
            "occupied_band_1based":list(occupied),"empty_band_1based":list(empty),
            "sampled_direct_gap_eV":float(eig[min(empty)-1]-eig[max(occupied)-1]),
            "pseudo_crossgap_transverse_spin_block_frobenius_norm":magnitude,
            "candidate_pair_count":len(pairs),"pairs":pairs}
    return {
      "status":"ZERO_JOB_PSEUDO_INTERBAND_SPIN_VERTEX_PROTOTYPE__NOT_PHYSICAL_KUBO",
      "input_sha256":{"WAVECAR":whash,"EIGENVAL":ehash},
      "input_names":{"WAVECAR":Path(wavecar).name,"EIGENVAL":Path(eigenval).name},
      "eigenval_header":meta,
      "source":"read-only raw pseudo spinors from already completed SOC VASP",
      "limitations":["PAW augmentation/operator overlap absent",
                     "pseudospin normalized in Euclidean pseudo metric for screening only",
                     "no velocity operator, position connection, or magnetic source derivative",
                     "no layer-odd spin projection, no whole-BZ integration",
                     "finite matrix norm is not a finite static NMEE"],
      "new_PBS_jobs":0,"material_static_NMEE":None,"k_points":results}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--wavecar",required=True)
    p.add_argument("--eigenval",required=True)
    p.add_argument("--expected-wavecar-sha256")
    p.add_argument("--expected-eigenval-sha256")
    p.add_argument("--out",required=True,help="new sanitized JSON path; refuses overwrite")
    a=p.parse_args()
    dest=Path(a.out)
    if dest.exists():raise FileExistsError("Cannot overwrite existing output: "+str(dest))
    report=screen(a.wavecar,a.eigenval,expected_wavecar_sha256=a.expected_wavecar_sha256,expected_eigenval_sha256=a.expected_eigenval_sha256)
    dest.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"k_points":list(report["k_points"]),
          "input_sha256":report["input_sha256"],"new_PBS_jobs":0},indent=2))
