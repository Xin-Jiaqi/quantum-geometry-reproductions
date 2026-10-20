"""Numerical multi-start global-gap screen on the reciprocal torus, not a proof.

Each sample has a grid search followed by local continuous minimization for
occupied maximum and empty minimum. Distinct initial valleys and multiple
coarse low-gap seeds test that the Dirac region is not hiding a Fermi pocket.
"""
import json, numpy as np
from scipy.optimize import minimize
from dataclasses import replace
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from scripts.run_stage3_sign_reversal import spec_at

def gap_check(u,t,n=46):
    spec=replace(spec_at(u),tunneling_eV=t)
    obj=SlidingBilayer(spec);G=obj.model().reciprocal
    f=np.linspace(0,1,n,endpoint=False)
    uv=np.stack(np.meshgrid(f,f,indexing='ij'),axis=-1).reshape(-1,2)
    ek=np.linalg.eigvalsh(obj.matrices(uv@G))
    distv=spec.mu_eV-ek[:,3];distc=ek[:,4]-spec.mu_eV
    outcomes=[]
    for key,sign in [('valence',-1),('conduction',+1)]:
        arr=distv if sign==-1 else distc
        starts=np.argsort(arr)[:12]
        initial=[uv[z] for z in starts]+[np.array([1/3,2/3]),np.array([2/3,1/3])]
        def val(q):
            en=np.linalg.eigvalsh(obj.matrices((np.asarray(q)%1)[None,:]@G))[0]
            return float((spec.mu_eV-en[3]) if sign==-1 else (en[4]-spec.mu_eV))
        loc=[]
        for seed in initial:
            opt=minimize(val,seed,method='Nelder-Mead',options={'xatol':2e-12,'fatol':2e-13,'maxiter':600})
            loc.append((float(opt.fun),list((opt.x%1).astype(float))))
        winner=min(loc,key=lambda x:x[0])
        outcomes.append({'band_edge':key,'minimum_distance_eV':winner[0], 'fractional_k':winner[1],
                         'grid_min_eV':float(arr.min()),'multiple_start_count':len(loc),
                         'max_negative_search_gap_eV':float(min([x[0] for x in loc]))})
    return {'u':u,'t_eV':t,'gap_to_mu_eV':min(x['minimum_distance_eV'] for x in outcomes),
            'full_gap_lower_bound_from_minima_eV':sum(x['minimum_distance_eV'] for x in outcomes),
            'edges':outcomes,'method':'grid + multi-start continuous Nelder-Mead local minimization; cannot mathematically certify global minimum'}

if __name__=='__main__':
    out={'cases':[gap_check(.2,t) for t in [.010,.015,.01963385,.020]]}
    p=ROOT/'results/sliding_stage5_global_gap_screen.json'
    p.write_text(json.dumps(out,indent=2)+'\n')
    for c in out['cases']:
        print('t',c['t_eV'],'mu min',c['gap_to_mu_eV'],'edge gap',c['full_gap_lower_bound_from_minima_eV'])
