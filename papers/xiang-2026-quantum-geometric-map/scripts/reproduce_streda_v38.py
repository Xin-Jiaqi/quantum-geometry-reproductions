#!/usr/bin/env python3
"""Independent clean two-band sea/metal boundary benchmark; not full PHE."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from qgmap.streda_clean_sea_v38 import bulk_sea,qwz_model,fhs_chern,metallic_pole_audit

def evaluate(mesh=61,mass=-1,mu=0,eta=.02):
    data=bulk_sea(mesh,mass,mu,eta)
    h,_,_=qwz_model(mesh,mass)
    result={'scientific_status':'PARTIAL',
      'scope':'clean two-band antisymmetric Kubo-Streda Fermi-sea benchmark ONLY',
      'units':'sigma/(e^2/h), hbar=e=1',
      'fhs_chern':fhs_chern(h),'sea':data,
      'invalid_inference':'not full metallic/disordered conductivity and not published EB PHE'}
    if abs(mu)>=.8:
        result['metallic_boundary_diagnostic']=metallic_pole_audit(mesh,mass,mu,eta)
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',type=int,default=61)
    ap.add_argument('--mass',type=float,default=-1.)
    ap.add_argument('--mu',type=float,default=0.)
    ap.add_argument('--eta',type=float,default=.02)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    x=json.dumps(evaluate(a.mesh,a.mass,a.mu,a.eta),indent=2,allow_nan=False)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(x+'\n')
    else:
        print(x)
if __name__=='__main__':main()
