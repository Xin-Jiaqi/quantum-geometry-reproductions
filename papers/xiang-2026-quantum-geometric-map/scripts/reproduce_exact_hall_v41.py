#!/usr/bin/env python3
"""v41 exact finite-torus linear Hall reference; original PRB EB remains OPEN."""
import argparse
import json
from pathlib import Path
from qgmap.exact_hall_v41 import audit

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--L',type=int,default=4)
    p.add_argument('--mass',type=float,default=-1.)
    p.add_argument('--W',type=float,default=.12)
    p.add_argument('--seed',type=int,default=1947)
    p.add_argument('--ntheta',type=int,default=7)
    p.add_argument('--mu',type=float,default=0.)
    p.add_argument('--eta',type=float,default=0.,help='response broadening, not SCBA eta')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    out=audit(a.L,a.mass,a.W,a.seed,a.ntheta,a.mu,a.eta)
    data=json.dumps(out,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(data)
    else:print(data,end='')
if __name__=='__main__':main()
