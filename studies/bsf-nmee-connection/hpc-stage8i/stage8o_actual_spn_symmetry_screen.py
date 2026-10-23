#!/usr/bin/env python3
"""Stage8o independent EIG/SPN symmetry diagnostics using actual VASP geometry.

READ-ONLY, zero PBS jobs. Requires numpy/spglib already available in harness.
Writes sanitized dimensionless checks only; never copies POTCAR/SPN/raw WAVECAR.
Validates actual crystallographic fractional Seitz rotations and ALL 441 k
points, unlike the Stage8n 12-pair wrong Cartesian mirror check.

Limits: 96-band SPN subspace may not be closed under physical symmetries.
This diagnostics cannot by itself establish *why* a covariance test fails;
sewing matrix and operator-complete Kubo checks remain separate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from stage8m_response import read_win_kpoints, read_win_cell, read_spn_formatted
from stage8m_plaq import read_eig


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def parse_poscar(path):
    lines=Path(path).read_text().splitlines()
    if len(lines)<14:raise ValueError("POSCAR truncated")
    scale=float(lines[1].split()[0])
    if scale<=0:raise ValueError("Only positive POSCAR scale supported")
    lattice=np.array([[float(s)*scale for s in row.split()[:3]]
                     for row in lines[2:5]])
    species=lines[5].split()
    counts=[int(v) for v in lines[6].split()]
    if species!=["W","Se"] or counts!=[2,4]:
        raise ValueError("Expected W2Se4 POSCAR (W, Se)")
    j=7
    if lines[j].strip().lower().startswith("s"):j+=1
    if not lines[j].strip().lower().startswith("d"):
        raise ValueError("Expected direct fractional positions")
    frac=np.array([[float(x) for x in row.split()[:3]]
                   for row in lines[j+1:j+7]])
    if frac.shape!=(6,3):raise ValueError("Bad six-atom site list")
    types=[1]*2+[2]*4
    return lattice,frac,types


def reciprocal_grid(kfrac,n1=21,n2=21,tol=2e-5):
    k=np.asarray(kfrac,float)
    if k.shape!=(n1*n2,2):raise ValueError("Expected complete two-dimensional uniform mesh")
    entries=np.rint(k*np.array([n1,n2])).astype(int)
    actual=(entries/np.array([n1,n2]))
    dif=np.abs(((k-actual+.5)%1)-.5)
    if float(np.max(dif))>tol:raise ValueError("k-grid does not match declared mesh")
    dic={}
    for ik,(i,j) in enumerate(entries):
        xy=(int(i%n1),int(j%n2))
        if xy in dic:raise ValueError("Duplicate folded fractional k")
        dic[xy]=ik
    if len(dic)!=n1*n2:raise ValueError("incomplete mesh")
    return dic


def pairs_from_fractional_symmetry(kfrac,R3,n1=21,n2=21):
    """Map ALL actual k to R^-T k modulo reciprocal integer lattice."""
    dic=reciprocal_grid(kfrac,n1,n2)
    R=np.asarray(R3,dtype=int)
    if R.shape!=(3,3) or not np.all(R[:2,2]==0) or not np.all(R[2,:2]==0):
        raise ValueError("Symmetry mixes k_z with in-plane k; requires a different mesh")
    RK=np.linalg.inv(R).T
    if not np.allclose(RK,np.rint(RK),atol=1e-12):
        raise ValueError("Non-integer symmetry reciprocal transformation")
    RK=np.rint(RK).astype(int)
    k=np.asarray(kfrac,float)
    pairs=[]
    for ik,(u,v) in enumerate(k):
        transformed=RK@np.array([u,v,0.])
        x=int(np.rint(transformed[0]*n1))%n1
        y=int(np.rint(transformed[1]*n2))%n2
        dest=dic[(x,y)]
        pairs.append((ik,dest))
    if len(set(j for _,j in pairs))!=n1*n2:
        raise AssertionError("Actual symmetry is not bijective on input mesh")
    return pairs


def spin_axial_Q(lattice,R3):
    """R3 is spglib fractional DIRECT lattice operation."""
    A=np.asarray(lattice,float).T
    rc=A@np.asarray(R3,float)@np.linalg.inv(A)
    if not np.allclose(rc.T@rc,np.eye(3),atol=1e-5):
        raise ValueError("Nonorthogonal Cartesian operation for actual cell")
    det=np.linalg.det(rc)
    return det*rc


def spin_gram_at_k(S):
    """S is [nb, nb, 3]. C_ij=Re Tr[S_i^dag S_j]."""
    return np.array([[np.vdot(S[:,:,i],S[:,:,j]).real for j in range(3)]
                     for i in range(3)])


def symmetry_summary(kfrac,energies,spins,lattice,R3,tfrac,n1=21,n2=21):
    pairs=pairs_from_fractional_symmetry(kfrac,R3,n1,n2)
    Q=spin_axial_Q(lattice,R3)
    E=np.asarray(energies,float)
    if E.ndim!=2 or E.shape[0]!=len(pairs):
        raise ValueError("Energy shape mismatch")
    S=np.asarray(spins)
    if S.shape[:1]!=(len(pairs),) or S.shape[-1]!=3:
        raise ValueError("SPN shape mismatch")
    C=[spin_gram_at_k(x) for x in S]
    normed=[]
    absres=[]
    spectrum=[]
    for i,j in pairs:
        ideal=Q@C[i]@Q.T
        err=np.linalg.norm(C[j]-ideal,"fro")
        denom=max(np.linalg.norm(C[j],"fro"),np.linalg.norm(ideal,"fro"),1e-12)
        normed.append(float(err/denom))
        absres.append(float(err))
        spectrum.append(float(np.max(np.abs(np.sort(E[i])-np.sort(E[j])))))
    R=np.asarray(R3,int)
    horizontal_glide=(np.array_equal(R,np.diag([1,1,-1])) and
                       np.any(np.abs(np.asarray(tfrac,float)[:2]%1)>1e-6))
    report={
        "direct_fractional_R":R.tolist(),
        "fractional_translation":np.asarray(tfrac,float).tolist(),
        "spin_axial_rotation_cartesian_Q":Q.tolist(),
        "mesh_pairs_exact":len(pairs),
        "max_bandwise_spectrum_symmetry_residual_eV":max(spectrum),
        "max_spin_gram_covariance_residual":max(normed),
        "median_spin_gram_covariance_residual":float(np.median(normed)),
        "max_spin_gram_covariance_absolute":max(absres),
        "max_global_spin_gram_norm":max(float(np.linalg.norm(Ci,"fro")) for Ci in C),
        "is_horizontal_glide":bool(horizontal_glide),
        "caveat":"Finite 96-band subspace may be non-closed at top; Gram covariance failure alone does NOT localize defect to VASP, ISYM=0 or MMN."
    }
    if horizontal_glide:
        # An actual g_z=mirror glide leaves in-plane k fixed, flips axial Sx,Sy.
        # Isolated eigenstates (gap to both neighbours) have zero spin x/y.
        good=[]
        maxspin=0.0
        skipped=0
        top=[]
        for ik in range(E.shape[0]):
            for ib in range(1,E.shape[1]-1):
                sep=min(abs(E[ik,ib]-E[ik,ib-1]),
                        abs(E[ik,ib+1]-E[ik,ib]))
                if sep < .001:
                    skipped+=1
                    continue
                sp=float(np.hypot(S[ik,ib,ib,0].real,S[ik,ib,ib,1].real))
                maxspin=max(maxspin,sp)
                if sp>1e-4:
                    good.append({"k_index":ik,"band_1based":ib+1,
                                 "local_band_gap_eV":float(sep),
                                 "Sxy_diag_norm_raw_Pauli":sp})
        report["glide_nondegenerate_transverse_spin_check"]={
            "threshold_band_isolation_eV":.001,
            "threshold_transverse_spin_raw_Pauli":1e-4,
            "degenerate_nearby_states_skipped":skipped,
            "max_transverse_diag_spin":maxspin,
            "number_of_above_threshold":len(good),
            "examples":sorted(good,key=lambda x:-x["Sxy_diag_norm_raw_Pauli"])[:12],
            "limitation":"Band96 has no 97 neighbour; excluded. Raw SPN diagonal is a pseudopotential-derived band-spin proxy, not a full spatial layer operator."
        }
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--poscar",required=True);p.add_argument("--win",required=True)
    p.add_argument("--eig",required=True);p.add_argument("--spn",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--symprec",type=float,default=.001)
    args=p.parse_args()
    out=Path(args.out)
    if out.exists():raise FileExistsError("Refuse to overwrite independent sanitized evidence")
    import spglib
    lattice,sites,types=parse_poscar(args.poscar)
    dataset=spglib.get_symmetry((lattice,sites,types),symprec=args.symprec)
    if not dataset:raise ValueError("spglib found no symmetries")
    kpts=read_win_kpoints(args.win)
    E=read_eig(args.eig)
    if E.shape[0]!=len(kpts):raise ValueError("EIG/WIN count discrepancy")
    S=read_spn_formatted(args.spn,E.shape[1],len(kpts))
    groups={}
    for R,t in zip(dataset["rotations"],dataset["translations"]):
        key=tuple(np.asarray(R).ravel())
        if key in groups:continue
        groups[key]=symmetry_summary(kpts,E,S,lattice,R,t)
    # Explicitly list SQUARE of differences of the 3 physical spin
    # components, including the correct axial spin rotation Q.
    payload={
      "status":"CORRECTED_EXACT_SYMMETRY_SCREEN__ROOT_CAUSE_NOT_ASSIGNED",
      "input_sha256":{name:sha256(getattr(args,name)) for name in ("poscar","win","eig","spn")},
      "spglib_symprec_A":args.symprec,
      "n_kpts":len(kpts),"n_bands":E.shape[1],
      "operations":list(groups.values()),
      "no_new_pbs_jobs":True,"NMEE_computed":False,
      "scientific_warning":"Corrected Gram tests alone cannot establish true glide sewing, full 96-band closure, MMN compatibility or trustworthy static response. No proof VASP ISYM=0 is responsible."
    }
    out.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps({"ops":len(groups),"done":True,"no_new_pbs_jobs":True},indent=2))


if __name__=="__main__":
    main()
