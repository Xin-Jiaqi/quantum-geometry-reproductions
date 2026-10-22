#!/usr/bin/env python3
"""Stage8l: zero-DFT pseudo-spin INTERBAND vertex pilot from existing WAVECAR.

Computes plane-wave pseudo <v|Sx/y|c> at selected k (no velocity, no PAW,
no layer operator, no E-field, no finite B, no NMEE). Purpose: falsifiable
pre-Wannier screen of occupied-empty spin-flip amplitudes. Requires numpy.
"""
import argparse
import json
from pathlib import Path
import numpy as np

from stage8k_independent_reader import (read_header,kpoint_meta,band_coeffs,
                                        split_spinor,find_kindex)
from stage8k_spin_analysis import parse_eigenval

# Fractional in actual two-dimensional reciprocal primitive basis
TARGETS={
 "GAMMA_M_TWO_THIRDS_NOT_KPRIME":(-1/3,1/3,0),
 "NEAR_M_NOT_EXACT_M":(10/21,0,0),
 "GAMMA_K_INTERIOR_NOT_CBM":(3/21,3/21,0),
 "TRUE_K":(1/3,1/3,0),
}
# Spinor band indices are explicit ONE-BASED here; NELECT=36.
OCCUPIED_BANDS=(33,34,35,36)
EMPTY_BANDS=tuple(range(37,45))

def source_vertices(pseudo_valence,pseudo_conduction):
    """Input (up,down) pseudo spinors. Return 1/2-Pauli offdiag vertices.

    Euclidean pseudo norm convention for DIAGNOSTIC ONLY.
    Do NOT equate to PAW physical spin or a Kubo Berry curvature.
    """
    vu,vd=pseudo_valence;cu,cd=pseudo_conduction
    norms=[np.vdot(x,x).real+np.vdot(y,y).real
           for x,y in (pseudo_valence,pseudo_conduction)]
    if min(norms)<=1e-10:raise ValueError("Invalid pseudo spinor norm")
    den=(norms[0]*norms[1])**.5
    ud=np.vdot(vu,cd);du=np.vdot(vd,cu)
    sx=(ud+du)/(2*den)
    sy=(-1j*ud+1j*du)/(2*den)
    sz=(np.vdot(vu,cu)-np.vdot(vd,cd))/(2*den)
    return np.array([sx,sy,sz])

def screen(wavecar,eigenval,occupied=OCCUPIED_BANDS,empty=EMPTY_BANDS):
    h=read_header(wavecar)
    kpts,energies=parse_eigenval(eigenval)
    if h["nk"]!=len(kpts) or h["nb"]!=len(energies[0]):
        raise ValueError("WAVECAR/EIGENVAL k/band counts not equal")
    if max(occupied)>=min(empty) or min(occupied)<1 or max(empty)>h["nb"]:
        raise ValueError("Occupied-empty window invalid")
    karr=np.asarray(kpts,float)
    results={}
    for name,coord in TARGETS.items():
        idx,distance=find_kindex(h,coord,karr)
        if distance>1e-5:
            raise ValueError(f"Requested k not represented exactly in k mesh: {name}, nearest distance {distance}")
        nplane,wk,enocc=kpoint_meta(wavecar,h,idx)
        eig=np.array([x[0] for x in energies[idx]])
        mismatch=float(np.max(np.abs(enocc[:,0]-eig)))
        if mismatch>1e-5:raise ValueError(f"Wavecar/Eigenval energy mismatch: {mismatch}")
        components={}
        bands=sorted(set(occupied).union(empty))
        for band in bands:
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
        magnitude=float(np.sqrt(sum(p["pseudo_abs_Sx"]**2+p["pseudo_abs_Sy"]**2
                                    for p in pairs)))
        results[name]={
            "requested_fractional_k":list(coord),
            "sampled_fractional_k":[float(a) for a in kpts[idx]],
            "k_index_0based":idx,"source_energy_roundtrip_max_diff_eV":mismatch,
            "occupied_band_1based":list(occupied),"empty_band_1based":list(empty),
            "sampled_direct_gap_eV":float(eig[min(empty)-1]-eig[max(occupied)-1]),
            "pseudo_crossgap_transverse_spin_block_frobenius_norm":magnitude,
            "candidate_pair_count":len(pairs),"pairs":pairs}
    return {
      "status":"ZERO_JOB_PSEUDO_INTERBAND_SPIN_VERTEX_PROTOTYPE__NOT_PHYSICAL_KUBO",
      "input_wavecar_name":Path(wavecar).name,
      "input_eigenval_name":Path(eigenval).name,
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
    p.add_argument("--out",required=True,help="new LOCAL sanitized JSON path; refuses overwrite")
    a=p.parse_args()
    dest=Path(a.out)
    if dest.exists():raise FileExistsError("Cannot overwrite existing output: "+str(dest))
    report=screen(a.wavecar,a.eigenval)
    dest.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],
          "k_points":list(report["k_points"]),"new_PBS_jobs":0},indent=2))
