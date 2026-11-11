#!/usr/bin/env python3
"""Stage8i: conservative analysis of a noncollinear spinor VASP EIGENVAL grid.

Sampled 21x21x1 BZ results are NOT proof of continuous-grid insulating gap.
Occupancy and source matrix-element tests are separate required gates.
"""
import argparse,json,math
from pathlib import Path

def analyze_eigenval(path,expected_kpoints=441):
    lines=Path(path).read_text(errors="replace").splitlines()
    if len(lines)<8:raise ValueError("EIGENVAL incomplete")
    meta=lines[5].split()
    if len(meta)<3:raise ValueError("EIGENVAL counts unavailable")
    nelect,nk,nb=(float(meta[0]),int(meta[1]),int(meta[2]))
    if nk!=expected_kpoints:
        raise RuntimeError(f"Wrong k mesh: {nk} != expected {expected_kpoints}")
    if not math.isfinite(nelect) or abs(nelect-round(nelect))>1e-7:
        return {"status":"UNVERIFIED_NONINTEGER_NELECT",
                "NELECT":nelect,"num_kpoints":nk,"num_bands":nb}
    occ=int(round(nelect))
    if occ<1 or occ>=nb or nb-occ<12:
        return {"status":"BANDS_INSUFFICIENT_OR_WRONG_SPINOR_COUNT",
                "NELECT":nelect,"nocc_spinor":occ,
                "num_kpoints":nk,"num_bands":nb}
    i=6
    valence_max=-float("inf")
    conduction_min=float("inf")
    min_direct=float("inf")
    valence_k=None
    conduction_k=None
    fractional_occupation_near_gap=False
    for k in range(nk):
        while i<len(lines) and not lines[i].strip():i+=1
        if i>=len(lines):raise ValueError("Missing EIGENVAL k-point")
        krow=lines[i].split()
        if len(krow)!=4:raise ValueError("EIGENVAL expected fractional k and weight")
        kfrac=[float(v) for v in krow[:3]]
        i+=1
        energies=[]
        occupations=[]
        for j in range(nb):
            if i>=len(lines):raise ValueError("Truncated bands")
            vals=lines[i].split();i+=1
            if len(vals)!=3:
                raise ValueError("Spinor EIGENVAL expected exactly band,energy,occupancy columns; not a collinear spin-two-band export")
            if int(vals[0])!=j+1:raise ValueError("Noncontiguous band indices")
            energies.append(float(vals[1]));occupations.append(float(vals[2]))
        v=energies[occ-1];c=energies[occ]
        if not all(math.isfinite(x) for x in (v,c)):
            raise ValueError("NaN/Inf band energies")
        if v>valence_max:valence_max=v;valence_k=kfrac
        if c<conduction_min:conduction_min=c;conduction_k=kfrac
        min_direct=min(min_direct,c-v)
        if occupations[occ-1]<.95 or occupations[occ]>.05:
            fractional_occupation_near_gap=True
    indirect=conduction_min-valence_max
    verdict=("SAMPLED_POSITIVE_GAP_WITH_SMEARING_CAVEAT" if indirect>0
             else "SAMPLED_OVERLAPPING_BANDS_OR_METALLIC")
    if fractional_occupation_near_gap:
        verdict="OCCUPATION_AMBIGUOUS_AT_0p05eV_SMEARING"
    return {"status":verdict,"NELECT":nelect,"nocc_spinor":occ,
        "num_kpoints":nk,"num_bands":nb,
        "sampled_indirect_gap_eV":indirect,
        "sampled_min_direct_gap_eV":min_direct,
        "valence_max_eV":valence_max,"conduction_min_eV":conduction_min,
        "VBM_fractional_k":valence_k,"CBM_fractional_k":conduction_k,
        "fractional_occupation_warning":fractional_occupation_near_gap,
        "assumptions":"noncollinear single-spinor EIGENVAL occupation max=1; spinor ordering consistent with actual NELECT; k-grid sampled only",
        "cannot_claim":["continuous 2D BZ gap","nonzero transverse interband spin mixing",
                        "finite NMEE response","Wannier window completeness"]}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("file");p.add_argument("--nk",type=int,default=441)
    a=p.parse_args()
    print(json.dumps(analyze_eigenval(a.file,a.nk),indent=2))
