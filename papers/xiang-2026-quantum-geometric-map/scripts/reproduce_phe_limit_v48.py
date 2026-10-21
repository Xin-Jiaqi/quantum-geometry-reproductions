#!/usr/bin/env python3
"""v48 symmetric spin-PHE: eta^2, mu/T and cutoff gates; author SM OPEN."""
import argparse,json
from pathlib import Path
from qgmap.phe_limit_audit_v48 import full_audit,eta_limit
from qgmap.dirac_mixed_rta_v46 import Spec

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['all','eta'],default='all')
    p.add_argument('--beta',type=float,default=.125)
    p.add_argument('--mu',type=float,default=.3)
    p.add_argument('--temperature',type=float,default=.03)
    p.add_argument('--nr',type=int,default=140)
    p.add_argument('--nphi',type=int,default=160)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    obj=full_audit() if args.mode=='all' else eta_limit(Spec(beta=args.beta,mu=args.mu,temperature=args.temperature,nr=args.nr,nphi=args.nphi))
    out=json.dumps(obj,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(out,encoding='utf-8')
    else:print(out,end='')
if __name__=='__main__':main()
