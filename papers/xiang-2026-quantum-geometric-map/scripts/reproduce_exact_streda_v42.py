#!/usr/bin/env python3
"""v42 finite-system Streda I/II and exact Lehmann: explicitly NOT original EB PHE."""
import argparse,json
from pathlib import Path
from qgmap.exact_streda_v42 import streda_audit,ensemble_audit

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--L',type=int,default=4)
    p.add_argument('--mass',type=float,default=-1.)
    p.add_argument('--W',type=float,default=.12)
    p.add_argument('--mu',type=float,default=1.1)
    p.add_argument('--seed',type=int,default=1947)
    p.add_argument('--theta-x',type=float,default=.11)
    p.add_argument('--theta-y',type=float,default=.33)
    p.add_argument('--etas',type=float,nargs='+',default=[.16,.05,.01,.003])
    p.add_argument('--ensemble-seeds',type=int,nargs='*',default=[])
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    th=(a.theta_x,a.theta_y)
    result=(ensemble_audit(a.L,a.mass,a.W,tuple(a.ensemble_seeds),th,a.mu,a.etas[-1])
       if a.ensemble_seeds else streda_audit(a.L,a.mass,a.W,a.seed,th,a.mu,tuple(a.etas)))
    txt=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(txt)
    else:print(txt,end='')
if __name__=='__main__':main()
