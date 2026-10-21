#!/usr/bin/env python3
"""v47 corrected original spin PHE (XY+YX)/2; paper SM matching remains OPEN."""
import argparse,json
from pathlib import Path
from qgmap.dirac_mixed_rta_v46 import Spec,decompose
from qgmap.dirac_source_tensor_v47 import tensor,derivative_source_check

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--beta',type=float,default=.125)
    p.add_argument('--mu',type=float,default=.3)
    p.add_argument('--temperature',type=float,default=.03)
    p.add_argument('--eta',type=float,default=.002)
    p.add_argument('--nr',type=int,default=160)
    p.add_argument('--nphi',type=int,default=192)
    p.add_argument('--cutoff',type=float,default=2.5)
    p.add_argument('--step',type=float,default=2e-5)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    spec=Spec(a.beta,a.mu,a.temperature,a.eta,a.cutoff,a.nr,a.nphi)
    v=tensor(spec,step=a.step)
    v['independent_full_Fermi_source_fd']=derivative_source_check(spec,step=a.step)
    out=json.dumps(v,allow_nan=False,indent=2)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(out,encoding='utf-8')
    else:print(out,end='')
if __name__=='__main__':main()
