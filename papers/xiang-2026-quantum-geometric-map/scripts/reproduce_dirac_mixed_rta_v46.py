#!/usr/bin/env python3
"""Two distinct mixed-field RTA collision targets; neither identified as paper SM."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.dirac_mixed_rta_v46 import Spec,audit,unprepared_two_field
from qgmap.dirac_phe_protocol_v45 import printed_eq11_geometry

def evaluate(spec=Spec(),h=3e-5,step=3e-5):
    out=audit(spec,h)
    other=unprepared_two_field(spec,step)
    out['unprepared_frozen_H0_collision_target']=other
    out['unprepared_tau_linear_coefficient']=spec.eta*other['total']
    out['printed_main_eq8_positive_geometry_reference']=-printed_eq11_geometry(spec.beta)*np.tanh(spec.mu/(2*spec.temperature))
    out['relative_difference_unprepared_vs_eq8']=out['unprepared_tau_linear_coefficient']/out['printed_main_eq8_positive_geometry_reference']-1
    out['important_note']='Formal Eq8 geometry is a separate reference; no author SM collision prescription matched.'
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name,default,tp in [('beta',.125,float),('mu',.3,float),
      ('temperature',.055,float),('eta',.065,float),('cutoff',2.5,float),
      ('nr',96,int),('nphi',144,int)]:
        ap.add_argument('--'+name,type=tp,default=default)
    ap.add_argument('--h',type=float,default=3e-5)
    ap.add_argument('--step',type=float,default=3e-5)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    spec=Spec(args.beta,args.mu,args.temperature,args.eta,args.cutoff,args.nr,args.nphi)
    txt=json.dumps(evaluate(spec,args.h,args.step),indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(txt)
    else:print(txt,end='')

if __name__=='__main__':main()
