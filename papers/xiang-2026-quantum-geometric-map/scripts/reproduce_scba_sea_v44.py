#!/usr/bin/env python3
"""Provisional SCBA Hall I(a)+II versus finite exact v42; NOT original EB response."""
import argparse,json
from pathlib import Path
from qgmap.scba_sea_v44 import compare_candidate,candidate_total

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--mesh',type=int,default=35)
    a.add_argument('--L',type=int,default=4)
    a.add_argument('--mass',type=float,default=-1.)
    a.add_argument('--W',type=float,default=.12)
    a.add_argument('--mu',type=float,default=1.1)
    a.add_argument('--eta',type=float,default=.12)
    a.add_argument('--first-seed',type=int,default=1947)
    a.add_argument('--seeds',type=int,default=24)
    a.add_argument('--ntheta',type=int,default=7)
    a.add_argument('--lower',type=float,default=-12.)
    a.add_argument('--quadrature',type=int,default=420)
    a.add_argument('--output',type=Path)
    args=a.parse_args()
    if args.seeds<=0:a.error('--seeds must be positive')
    kw=dict(mesh=args.mesh,L=args.L,mass=args.mass,W=args.W,
            mu=args.mu,eta=args.eta,lower=args.lower,quadrature=args.quadrature)
    result=(candidate_total(args.mesh,args.mass,0,args.mu,args.eta,lower=args.lower,
                           quadrature=args.quadrature) if args.W==0 else
             compare_candidate(**kw,seeds=tuple(range(args.first_seed,args.first_seed+args.seeds)),
                               ntheta=args.ntheta))
    out=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(out)
    else:print(out,end='')
if __name__=='__main__':main()
