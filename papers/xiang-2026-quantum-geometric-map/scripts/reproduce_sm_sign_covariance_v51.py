#!/usr/bin/env python3
"""v51 authenticated source: Eq20/31 relative sign under true spin-basis covariance."""
import argparse,hashlib,json
from pathlib import Path
from qgmap.sm_sign_covariance_v51 import analyze,integrated_table
SM_SHA='9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'
MAIN_SHA='b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--beta',type=float,default=.125)
    ap.add_argument('--mu',type=float,default=.3)
    ap.add_argument('--temperature',type=float,default=.03)
    ap.add_argument('--radial',type=int,default=80)
    ap.add_argument('--angular',type=int,default=112)
    ap.add_argument('--sm-pdf',type=Path)
    ap.add_argument('--main-pdf',type=Path)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    for filepath,expected in ((a.sm_pdf,SM_SHA),(a.main_pdf,MAIN_SHA)):
        if filepath and hashlib.sha256(filepath.read_bytes()).hexdigest()!=expected:
            ap.error('Original source SHA256 mismatch: '+str(filepath))
    if (a.beta,a.mu,a.temperature,a.radial,a.angular)==(.125,.3,.03,80,112):
        obj=analyze()
    else:obj=integrated_table(beta=a.beta,mu=a.mu,T=a.temperature,nr=a.radial,nphi=a.angular)
    obj['source_files']={'main_sha256':MAIN_SHA,'sm_sha256':SM_SHA,
          'main_bytes_checked':bool(a.main_pdf),'sm_bytes_checked':bool(a.sm_pdf),
          'main_source':'PRB Eq10','SM_source':'SM p5 Eq20; p8 Eqs31-33; p13 Eq62'}
    data=json.dumps(obj,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(data,encoding='utf8')
    else:print(data,end='')
if __name__=='__main__':main()
