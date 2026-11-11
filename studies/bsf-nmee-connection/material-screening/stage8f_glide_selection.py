"""Stage8f: exact symmetry of BiDB-2685 WSe2 AB-to-BA 3 shortest slide paths.

SOURCE: six-atom BiDB row 2685 z-scan, NOT fully ionic-relaxed.
This validates spatial selection rules for insulating *equilibrium/static*
Gamma_{i,z;j}, not a driven-DC NMEE nor a DFT response magnitude.

The total spin is layer-even; layer-difference spin is layer-odd.
A layer-exchange glide reverses total response, but not staggered response.
"""
import importlib.util
from pathlib import Path
import numpy as np

DIR=Path(__file__).resolve().parent
ROOT=DIR.parent
SPEC=importlib.util.spec_from_file_location("independent_motif",DIR/"ideal_stack_symmetry.py")
_motif=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(_motif)

ARMS={"A":(1/3,-1/3),"B":(-2/3,-1/3),"C":(1/3,2/3)}
TOP=(1,4,5)
LABELS=("W","W","Se","Se","Se","Se")
SOURCE=DIR/"source-structures"/"WSe2_BiDB_row2685_3R_polar_zscan.vasp"


def read_source():
    lines=SOURCE.read_text().splitlines()
    assert "row 2685" in lines[0] and "NOT full relaxation" in lines[0]
    cell=np.array([list(map(float,x.split())) for x in lines[2:5]])
    frac=np.array([list(map(float,x.split()[:3])) for x in lines[8:14]])
    frac[:,2] *= cell[2,2]
    frac[:,2] -= (frac[:,2].min()+frac[:,2].max())/2
    return cell,frac


def slide(frac,arm,s):
    if arm not in ARMS or not 0<=s<=1:
        raise ValueError("Bad arm or displacement fraction")
    p=frac.copy()
    p[list(TOP),:2]=(p[list(TOP),:2]+s*np.asarray(ARMS[arm]))%1
    return p


def spin_tensor_rank(ops,layer_odd=False):
    """Reynolds rank on 2x2 spin(i)-Zeeman(z)-E(j) static tensor.

    Under spatial z flip, total spin response picks f_z.
    Layer-odd spin adds another f_z, canceling it.
    """
    bases=[]
    for i in range(2):
        for j in range(2):
            seed=np.zeros((2,2));seed[i,j]=1
            projector=np.zeros((2,2))
            for _,_,fz,R2 in ops:
                projector += (1 if layer_odd else fz)*(R2@seed@R2.T)
            bases.append((projector/len(ops)).ravel())
    return int(np.linalg.matrix_rank(np.stack(bases),tol=1e-8))


def mirror_maps(left,right,cell):
    """Horizontal reflection between full six-atom motifs plus one xy shift."""
    return _motif.maps(left,list(LABELS),cell[:2,:2],
                       right,list(LABELS),np.eye(2),-1)


def gate():
    cell,src=read_source()
    cart=np.asarray(list(ARMS.values()))@cell[:2,:2]
    lengths=np.linalg.norm(cart,axis=1)
    results={}
    for arm in ARMS:
        rows=[]
        for s in (0,.25,.5,.75,1):
            f=slide(src,arm,s)
            ops=_motif.operations(f,list(LABELS),cell[:2,:2])
            rows.append({"s":s,"order":len(ops),
                         "total_spin_tensor_dim":spin_tensor_rank(ops),
                         "layer_odd_tensor_dim":spin_tensor_rank(ops,layer_odd=True),
                         "is_layer_exchange_glide":mirror_maps(f,f,cell),
                         "mirror_maps_to_one_minus_s":mirror_maps(f,slide(src,arm,1-s),cell)})
        results[arm]=rows
    return {"source":"BiDB row 2685 6-atom z-scan, NOT full ionic relaxation",
            "three_C3_related_arm_lengths_A":lengths.tolist(),
            "tensor":"M_ij=Gamma_{i,z;j}; equilibrium STATIC insulator",
            "paths":results,
            "prediction":"M_plus(s)=-M_plus(1-s); M_minus(s)=+M_minus(1-s) under exact layer exchange",
            "not_calculated":"real-material NMEE, DFT total energies, relaxed barrier, NEB, SOC full-BZ gap"}


if __name__=="__main__":
    import json
    print(json.dumps(gate(),indent=2))
