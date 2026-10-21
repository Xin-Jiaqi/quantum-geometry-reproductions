#!/usr/bin/env python3
"""Original SM source Eq4-10 finite-frequency and independently regularized EB pole."""
import argparse,hashlib,json
from pathlib import Path
from qgmap.sm_literal_spin_phe_v49 import SourceConfig
from qgmap.sm_two_color_liouville_v50 import Tone,two_color,run_limit

MAIN_SHA='b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730'
SM_SHA='9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4'
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name,value,tp in [('beta',.125,float),('mu',.3,float),('temperature',.03,float),
                          ('nr',100,int),('nphi',128,int),('cutoff',2.5,float),
                          ('eta-a',.002,float),('eta-b',.001,float),
                          ('omega-a',0.,float),('omega-b',0.,float),('dk',2e-5,float)]:
        p.add_argument('--'+name,type=tp,default=value)
    p.add_argument('--helicity',choices=['main','sm'],default='main')
    p.add_argument('--component',choices=['xy','yx'],default='xy')
    p.add_argument('--scan',action='store_true')
    p.add_argument('--sm-pdf',type=Path)
    p.add_argument('--main-pdf',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    for fn,sha in ((a.sm_pdf,SM_SHA),(a.main_pdf,MAIN_SHA)):
        if fn is not None and hashlib.sha256(fn.read_bytes()).hexdigest()!=sha:
            p.error(f'original source PDF SHA256 mismatch: {fn}')
    c=SourceConfig(a.beta,a.mu,a.temperature,a.nr,a.nphi,a.cutoff)
    t=Tone(a.omega_a,a.omega_b,a.eta_a,a.eta_b,a.dk)
    data=(run_limit(c,eta_b=a.eta_b,dk=a.dk,helicity=a.helicity,direction=a.component)
          if a.scan else two_color(c,t,a.helicity,a.component))
    data['source_provenance']={
        'main_pdf_sha256':MAIN_SHA,'supplemental_pdf_sha256':SM_SHA,
        'main_pdf_byte_verified':a.main_pdf is not None,'sm_pdf_byte_verified':a.sm_pdf is not None,
        'source_pages':'SM p3-5 Eqs4-20, p8 Eqs31-33; main p4 Eq10-11'}
    def serialize(v):
        if isinstance(v,complex):return {'real':v.real,'imag':v.imag}
        raise TypeError(type(v).__name__)
    out=json.dumps(data,default=serialize,indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(out,encoding='utf-8')
    else:print(out,end='')
if __name__=='__main__':main()
