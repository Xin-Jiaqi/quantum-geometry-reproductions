#!/usr/bin/env python3
"""Compare original SM literal Eq20 vs Eq31/32 and both printed Dirac helicities."""
import argparse, hashlib, json
from pathlib import Path
from dataclasses import asdict
from qgmap.sm_literal_spin_phe_v49 import SourceConfig,compare_source

SM_SHA='9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'
MAIN_SHA='b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--beta',type=float,default=.125)
    p.add_argument('--mu',type=float,default=.3)
    p.add_argument('--temperature',type=float,default=.03)
    p.add_argument('--nr',type=int,default=100)
    p.add_argument('--nphi',type=int,default=128)
    p.add_argument('--cutoff',type=float,default=2.5)
    p.add_argument('--main-pdf',type=Path,help='Optionally SHA256-verify the authentic main PDF')
    p.add_argument('--sm-pdf',type=Path,help='Optionally SHA256-verify the authentic SM PDF')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    for supplied,expected in ((a.main_pdf,MAIN_SHA),(a.sm_pdf,SM_SHA)):
        if supplied is not None and hashlib.sha256(supplied.read_bytes()).hexdigest()!=expected:
            p.error(f'uploaded original PDF SHA256 mismatch: {supplied}')
    cfg=SourceConfig(a.beta,a.mu,a.temperature,a.nr,a.nphi,a.cutoff)
    result=compare_source(cfg)
    result['config']=asdict(cfg)
    result['original_source_sha256']={'main_pdf':MAIN_SHA,'supplemental_pdf':SM_SHA}
    result['source_pdf_runtime_verified']={'main':a.main_pdf is not None,'sm':a.sm_pdf is not None}
    out=json.dumps(result,allow_nan=False,indent=2)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(out,encoding='utf-8')
    else:print(out,end='')
if __name__=='__main__':main()
