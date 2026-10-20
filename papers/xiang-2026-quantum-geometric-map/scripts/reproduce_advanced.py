#!/usr/bin/env python3
"""Reproducible advanced model entry point (explicitly *not* a real-material solver)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from qgmap.spin_multiband_v5 import four_band_audit
from qgmap.orbital_phe_v4 import orbital_phe_eq56
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d

CASES={
    "spin-four-band":(four_band_audit,23,"nondegenerate tau^0 density-matrix spin-MNHE, SM Eqs.16-19 vs 29"),
    "orbital-phe":(orbital_phe_eq56,61,"two-band 2D Bz PHE; literal SM Eq.56 vs field derivative"),
    "orbital-mnhe-2d":(orbital_mnhe_2d,61,"two-band strict 2D Bz squared Berry curvature only"),
    "orbital-mnhe-3d":(orbital_mnhe_3d,17,"3D two-band positional-shift vs literal SM Eq.52"),
}
LIMITS={
    "spin-four-band":"nondegenerate synthetic four-band; no impurity or material Wannier operator",
    "orbital-phe":"two-band periodic 2D; SM Eq.56 line disagreement OPEN; no omitted three-band terms",
    "orbital-mnhe-2d":"strict 2D Bz SBC sector; this does not test QMQ",
    "orbital-mnhe-3d":"two-band 3D; SM Eq.52 factor-two disagreement OPEN; not general multiband",
}

def evaluate(case,mesh=None):
    if case not in CASES: raise ValueError("Unknown case; use one of "+",".join(CASES))
    fn,n,desc=CASES[case]
    if mesh is not None:
        if mesh<9:raise ValueError("Mesh must be >=9")
        n=mesh
    result=fn(N=n)
    return {"schema":"qgmap-advanced-1","case":case,
            "science_status":"PARTIAL","validation_scope":"synthetic model / method subset",
            "description":desc,"limitations":LIMITS[case],"source":"v30 offline independently audited implementation",
            "parameters":{"N":n},"result":result}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",required=True,choices=sorted(CASES))
    parser.add_argument("--mesh",type=int,help="override recommended diagnostic mesh")
    parser.add_argument("--output",type=Path,help="write JSON instead of stdout")
    opt=parser.parse_args()
    data=json.dumps(evaluate(opt.case,opt.mesh),indent=2,allow_nan=False)
    if opt.output is None:print(data)
    else:
        opt.output.parent.mkdir(parents=True,exist_ok=True)
        opt.output.write_text(data+"\n",encoding="utf-8")
if __name__=="__main__":main()
