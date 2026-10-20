"""Static finite-field NMEE integral, independent of nondegenerate band sums.

This is an EQUILIBRIUM STATIC diagnostic, not the paper's driven DC branch.
"""
from __future__ import annotations
import argparse,sys,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from nmee_lab.clean_core.core import mesh,chi_k,fermi,MU_B


def response(spec:SlidingSpec, n=41,step=.25,offset=(.37,.23),batch=2048):
    obj=SlidingBilayer(spec);model=obj.model()
    k=mesh(model.reciprocal,n,offset=offset)
    area=abs(np.linalg.det(model.reciprocal))/(2*np.pi)**2/len(k)
    coeff={key:0. for key in ('top_x','bottom_x','total_x')}
    min_gap=1e12;occgap=1e12
    for i in range(0,len(k),batch):
        kr=k[i:i+batch];H0=obj.matrices(kr);V=obj.matrices(kr,axis=0)
        e0=np.linalg.eigvalsh(H0); f0=fermi(e0,spec.mu_eV,0.)
        min_gap=min(min_gap,float((e0[:,4]-e0[:,3]).min()))
        occgap=min(occgap,float(np.min(e0[:,4]-spec.mu_eV)),float(np.min(spec.mu_eV-e0[:,3])))
        val={key:[] for key in coeff}
        for sign in (+1,-1):
            B=sign*step
            h=H0-model.g_factor*MU_B*B*obj.spins['total_z']
            e,U=np.linalg.eigh(h)
            if not np.all((e[:,3]<spec.mu_eV)&(e[:,4]>spec.mu_eV)):
                raise ValueError('magnetic field closed Fermi gap')
            for key in coeff:
                v=chi_k(e,U,V,np.broadcast_to(obj.spins[key],H0.shape),f0)
                val[key].append(float(np.sum(v)))
        for key in coeff:
            coeff[key]+=area*(val[key][0]-val[key][1])/(2*step)
    return {'shift':list(spec.shift),'n':n,'B_step_T':step,'mu_eV':spec.mu_eV,
            'tunneling_eV':spec.tunneling_eV,'decay_nm':spec.decay_nm,
            'cross_sublattice':spec.cross_sublattice,'response':coeff,
            'sampled_middle_band_gap_eV':min_gap,'min_mu_to_band_eV':occgap,
            'total_vs_layers_error':coeff['total_x']-coeff['top_x']-coeff['bottom_x']}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,default=31)
    p.add_argument('--t',type=float,default=.028)
    p.add_argument('--cross-only',action='store_true')
    args=p.parse_args()
    s=SlidingSpec(tunneling_eV=args.t,cross_sublattice=not args.cross_only)
    start=time.time();print(json.dumps(response(s,n=args.n),indent=2));print('seconds',time.time()-start)

def response_weighted(spec, k, w, step=.25, batch=2048):
    """Arbitrary positive 2D quadrature with physical area weights."""
    obj=SlidingBilayer(spec);model=obj.model()
    k=np.asarray(k,float);w=np.asarray(w,float)
    if k.shape!=(len(w),2) or not np.all(np.isfinite(w)) or (w<0).any():
        raise ValueError('wrong weighted integration inputs')
    coeff={key:0. for key in ('top_x','bottom_x','total_x')}
    min_gap=1e9;occgap=1e9
    for i in range(0,len(k),batch):
        kr=k[i:i+batch];weight=w[i:i+batch];H0=obj.matrices(kr);V=obj.matrices(kr,axis=0)
        e0=np.linalg.eigvalsh(H0); f0=fermi(e0,spec.mu_eV,0.)
        min_gap=min(min_gap,float(np.min(e0[:,4]-e0[:,3])))
        occgap=min(occgap,float(np.min(e0[:,4]-spec.mu_eV)),float(np.min(spec.mu_eV-e0[:,3])))
        val={key:[] for key in coeff}
        for sign in (+1,-1):
            h=H0-sign*model.g_factor*MU_B*step*obj.spins['total_z']
            e,U=np.linalg.eigh(h)
            if not np.all((e[:,3]<spec.mu_eV)&(e[:,4]>spec.mu_eV)):
                raise ValueError('magnetic field closed gap')
            for key in coeff:
                raw=chi_k(e,U,V,np.broadcast_to(obj.spins[key],H0.shape),f0)
                val[key].append(float(np.dot(weight,raw)))
        for key in coeff:
            coeff[key]+=(val[key][0]-val[key][1])/(2*step)
    return {'response':coeff,'sampled_middle_band_gap_eV':min_gap,'min_mu_to_band_eV':occgap,
            'total_vs_layers_error':coeff['total_x']-coeff['top_x']-coeff['bottom_x']}
