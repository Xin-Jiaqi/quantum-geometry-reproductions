#!/usr/bin/env python3
"""Reproduce printed PRB Fig2(b) Eq11 conductance and ~43 microV Hall estimate."""
from pathlib import Path
import argparse,hashlib,json
from qgmap.fig2_si_v52 import Fig2Config,full_report,SOURCE_MAIN_SHA256,SOURCE_SM_SHA256

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--beta',type=float,default=.125)
    p.add_argument('--tau-ps',type=float,default=.5)
    p.add_argument('--g',type=float,default=10)
    p.add_argument('--bx-tesla',type=float,default=1)
    p.add_argument('--ey-vm',type=float,default=1e4)
    p.add_argument('--width-um',type=float,default=100)
    p.add_argument('--resistance-ohm',type=float,default=1000)
    p.add_argument('--main-pdf',type=Path,help='Optional original 8-page PRB PDF SHA256 check')
    p.add_argument('--sm-pdf',type=Path,help='Optional original 16-page SM SHA256 check')
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    for fn,expected in ((args.main_pdf,SOURCE_MAIN_SHA256),(args.sm_pdf,SOURCE_SM_SHA256)):
        if fn is not None and hashlib.sha256(fn.read_bytes()).hexdigest()!=expected:
            p.error('original PDF sha256 mismatch: '+str(fn))
    c=Fig2Config(args.beta,args.tau_ps,args.g,args.bx_tesla,args.ey_vm,args.width_um,args.resistance_ohm)
    result=full_report(c)
    result['source_byte_verified']={'main':args.main_pdf is not None,'sm':args.sm_pdf is not None}
    txt=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(txt,encoding='utf8')
    else:print(txt,end='')
if __name__=='__main__':main()
