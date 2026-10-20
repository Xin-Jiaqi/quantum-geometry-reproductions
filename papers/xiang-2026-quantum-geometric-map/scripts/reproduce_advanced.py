#!/usr/bin/env python3
"""Run named, protocol-separated model reproductions; not a real-material solver."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from qgmap.spin_multiband_v5 import four_band_audit
from qgmap.orbital_phe_v4 import orbital_phe_eq56
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d
from qgmap.warped_disorder_v17 import Lattice,ra_observable,density_ward,slope as ra_slope
from qgmap.spectral_kubo_v18 import spectral_bubble,slope as spectral_slope
from qgmap.finite_qomega_ward_v19 import finite_ward


def _spin(n): return four_band_audit(N=n)
def _phe(n): return orbital_phe_eq56(N=n)
def _mnhe2d(n): return orbital_mnhe_2d(N=n)
def _mnhe3d(n): return orbital_mnhe_3d(N=n)

def _scba_ra(n):
    p=Lattice(n=n)
    r=ra_observable(p,b=0.,protocol="onsite",vertex=True)
    onsite=ra_slope(p,step=.02,protocol="onsite",vertex=True)
    shift=ra_slope(p,step=.02,protocol="shift",vertex=True)
    ward=density_ward(p,b=0.,protocol="onsite")
    return {"protocol":"SCBA_RA_finite_eta_periodic_lattice",
            "model":"independent Wilson/warped periodic two-band regulator",
            "units":"dimensionless BZ mean; b is model energy or shift, not Tesla",
            "at_zero":{k:r[k] for k in ("symmetric","xy","yx","vertex_residual","scba_residual","scba_iterations")},
            "onsite_zeeman_slope":onsite,
            "shift_control_slope":shift,
            "q0_density_ward":{k:ward[k] for k in ("abs_error","residual","iterations","causality_max_eigenvalue")},
            "limit":"RA Fermi-surface bubble only; not full Kubo-Streda or published PHE"}

def _spectral(n):
    p=Lattice(n=n)
    r=spectral_bubble(p,b=0.,protocol="onsite",vertex=True)
    s=spectral_slope(p,step=.02,protocol="onsite",vertex=True)
    return {"protocol":"SCBA_symmetric_spectral_Kubo_four_channels",
            "model":"independent periodic two-band scalar-disorder regulator",
            "units":"dimensionless BZ mean; common e^2*hbar/pi and k measure omitted",
            "at_zero":r,"onsite_zeeman_slope":s,
            "limit":"RA+AR-RR-AA with matching noncrossing ladders; not full Kubo-Streda"}

def _ward(n):
    p=Lattice(n=n)
    r=finite_ward(p,q=.08,omega=.04,direction=1,b=0.,protocol="onsite")
    return {"protocol":r["protocol"],"model":"independent periodic two-band scalar-disorder regulator",
            "parameters":{k:r[k] for k in ("q","omega","eta","W","direction","field_protocol","b")},
            "residuals":r["residuals"],
            "vertex_norms":{k:r[k] for k in ("density_ladder_correction_norm","longitudinal_ladder_correction_norm","density_vertex_norm")},
            "iterations":r["iterations"],"warnings":r["warnings"],
            "limit":"finite-eta resolvent Ward identity; not physical continuity at eta>0"}

CASES={
    "spin-four-band":(_spin,23,"Four-band spin-MNHE density matrix versus SM Eq.29","Synthetic nondegenerate four-band model; no impurity or material Wannier operator."),
    "orbital-phe":(_phe,61,"2D Bz orbital PHE; literal SM Eq.56 versus finite-field derivative","Two-band periodic 2D model; printed Eq.56 discrepancy remains open."),
    "orbital-mnhe-2d":(_mnhe2d,61,"Strict 2D Bz orbital MNHE squared-curvature sector","Strict 2D SBC sector only; this does not test general QMQ."),
    "orbital-mnhe-3d":(_mnhe3d,17,"3D two-band positional shift versus literal SM Eq.52","Two-band 3D model; Eq.52 factor-two discrepancy remains open."),
    "scba-ra":(_scba_ra,41,"SCBA self-energy, RA bubble and noncrossing current ladder","Independent periodic regulator; RA Fermi-surface kernel only."),
    "kubo-spectral":(_spectral,31,"Symmetric four-channel spectral Kubo with RA, AR, RR, AA ladders","Zero-temperature finite-eta spectral sector; not complete Kubo-Streda."),
    "finite-qomega-ward":(_ward,31,"Finite-(q,omega) SCBA resolvent Ward identity","External eta is a regulator; real-frequency physical continuity and dc limits remain open."),
}

def evaluate(case,mesh=None):
    if case not in CASES: raise ValueError("Unknown case; use one of "+",".join(CASES))
    fn,default_mesh,description,limits=CASES[case]
    n=default_mesh if mesh is None else mesh
    if n<12: raise ValueError("Mesh must be >=12 for the shared CLI")
    result=fn(n)
    return {"schema":"qgmap-advanced-2","case":case,
            "science_status":"PARTIAL","validation_scope":"named synthetic model / method subset",
            "description":description,"limitations":limits,
            "source":"v30 offline independently audited implementation",
            "parameters":{"N":n},"result":result}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",required=True,choices=sorted(CASES))
    parser.add_argument("--mesh",type=int,help="override the case's recommended diagnostic mesh")
    parser.add_argument("--output",type=Path,help="write JSON instead of stdout")
    opt=parser.parse_args()
    data=json.dumps(evaluate(opt.case,opt.mesh),indent=2,allow_nan=False)
    if opt.output is None:print(data)
    else:
        opt.output.parent.mkdir(parents=True,exist_ok=True)
        opt.output.write_text(data+"\n",encoding="utf-8")
if __name__=="__main__":main()
