#!/usr/bin/env python3
"""Original PRB symmetric spin-PHE Eq(8)/Eq(11) vs full-plane RTA field shift."""
import argparse
import json
from pathlib import Path
from qgmap.dirac_phe_protocol_v45 import audit

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--beta',type=float,default=1/8)
    p.add_argument('--mu',type=float,default=.30)
    p.add_argument('--temperature',type=float,default=.055)
    p.add_argument('--zeeman-shift',type=float,default=.10)
    p.add_argument('--grid',type=int,default=301)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    if args.mu==0:p.error('mu=0 requires separate Dirac point protocol')
    out=audit(args.beta,args.mu,args.temperature,args.zeeman_shift,args.grid)
    content=json.dumps(out,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(content,encoding='utf8')
    else:print(content,end='')
if __name__=='__main__':main()
