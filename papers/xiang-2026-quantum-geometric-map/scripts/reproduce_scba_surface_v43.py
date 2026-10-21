#!/usr/bin/env python3
"""v43 reproducible SCBA/finite-exact Hall I(a) comparison, NOT full Hall."""
import argparse,json
from pathlib import Path
from qgmap.scba_surface_v43 import comparison

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--L',type=int,default=4)
    ap.add_argument('--mesh',type=int,default=35)
    ap.add_argument('--mass',type=float,default=-1.)
    ap.add_argument('--W',type=float,default=.12)
    ap.add_argument('--mu',type=float,default=1.1)
    ap.add_argument('--eta',type=float,default=.12)
    ap.add_argument('--first-seed',type=int,default=1947)
    ap.add_argument('--seeds',type=int,default=24)
    ap.add_argument('--ntheta',type=int,default=7)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    if a.seeds<1:ap.error('--seeds must be >=1')
    v=comparison(a.L,a.mesh,a.mass,a.W,a.mu,a.eta,
                 tuple(range(a.first_seed,a.first_seed+a.seeds)),a.ntheta)
    text=json.dumps(v,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text)
    else:print(text,end='')
if __name__=='__main__':main()
