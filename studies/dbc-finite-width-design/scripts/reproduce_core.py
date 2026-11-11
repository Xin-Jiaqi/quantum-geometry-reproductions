#!/usr/bin/env python3
"""Run calculations, never report a frozen JSON as a fresh evaluation."""
from __future__ import annotations
import argparse,json,sys,platform,time
from pathlib import Path
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from optical_design import FAMILIES,critical,predict,integrate_point,energy
from optical_design.models import h_only,kernels
from optical_design.lineshapes import optimum
from optical_design.two_feature import predict_two_feature
from optical_design.projectors import projector_kernel
from optical_design.closure import NormalForm,skew_optimum


def calculate(extended=False):
    shapes={s:optimum(s) for s in ('gaussian','lorentzian')}
    cases=[]
    for name in ('constant','bond_y','mixed_orbital'):
        p=FAMILIES[name];c=critical(p);pr=predict(c,1e-4)
        pp=replace(p,mass=pr['predicted_control'])
        ph=float(energy(0.,np.pi,pp)[0])+pr['photon_offset_from_saddle']
        full=integrate_point(pp,ph,1e-4)
        cases.append(dict(family=name,critical=c,prediction=pr,full_bz=full))
    local_errors=[]
    for name in ('constant','bond_y','mixed_orbital','dark'):
        p=FAMILIES[name];k=np.array([.31,2.23])
        for axes in ((0,0,1),(1,1,1),(0,0,0)):
            r=projector_kernel(lambda q:h_only(*q,p),k,axes=axes)
            expected=kernels(*k,p,axes)[1]
            local_errors.append(float(np.max(abs(expected-r['kernels']))))
    # Fresh independent full-BZ check of the previously discovered competitor.
    competitor=predict_two_feature(FAMILIES['mixed_orbital'],5e-5)
    pc=replace(FAMILIES['mixed_orbital'],mass=competitor['parameters']['mass'])
    competitor['fresh_full_bz']=integrate_point(pc,competitor['prediction']['photon'],5e-5)
    results=dict(status='MODEL_PROTOCOL_VALIDATED',paper='10.1103/sqwy-xhdf',
        source_type='Independent research extension, not author results',
        protocol='2D insulating fixed-basis three-orbital shift current; optical spectral convolution',
        units='Model energy in eV, dimensionless momentum; SI prefactor and thickness omitted',
        shapes=shapes,cases=cases,competitor=competitor,
        projector_max_absolute_error=max(local_errors),
        scope=['No real-material or DFT results','Not a microscopic collision model',
               'No generic crossing regularization or universal/global optimum theorem',
               'Photon frequency and control must both be optimized'],
        environment=dict(python=platform.python_version(),numpy=np.__version__))
    if extended:
        results['skew_scan']=[skew_optimum(nu,s) for s in ('gaussian','lorentzian') for nu in (0.,.5,1.,2.)]
        results['quadrature_check']=[]
        for case in cases:
            p=replace(FAMILIES[case['family']],mass=case['prediction']['predicted_control'])
            ref=integrate_point(p,case['full_bz']['photon'],1e-4,yorder=12,xorder=18,ypanels=48,xpanels=96)
            results['quadrature_check'].append(dict(family=case['family'],value=ref['total'],
                relative_change=ref['total']/case['full_bz']['total']-1))
    return results


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--extended',action='store_true')
    ap.add_argument('--output',type=Path)
    a=ap.parse_args();r=calculate(a.extended);text=json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text,encoding='utf-8')
    print(text,end='')

if __name__=='__main__':main()
