#!/usr/bin/env python3
"""Recompute model evidence; stdout contains JSON only. No author-curve fitting."""
import json
from dataclasses import replace
import numpy as np
from orbital_igmc import CuMnAs, geometry, integrate, mesh, response_si, density_reference


def solve(model,n=161,mu=0,T=25):
    geom=geometry(*model.operators(mesh(n)),group_sizes=model.group_sizes)
    return geom,{p:integrate(geom,mu,T,protocol=p) for p in ('printed_D10','density_A9')}


def serial(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def main():
    model=CuMnAs(a_angstrom=3.82)
    g,results=solve(model)
    benchmark={}
    for p,r in results.items():
        si=response_si(r,model,.0001,convention='harmonic')
        benchmark[p]={'integrated_kernel':r,
                      'imag_chi_harmonic_A_per_m_T':{k:v.imag for k,v in si.items()}}
    convergence=[]
    for n in (61,101,161,241):
        _,rs=solve(model,n)
        convergence.append({'mesh':n,'orbital_x_printed':rs['printed_D10']['orbital'][0],
                             'spin_x':rs['printed_D10']['spin'][0]})
    q=mesh(41);_,fast=solve(model,41)
    ref=density_reference(model,q,0,25)
    reference_error=max(float(np.max(np.abs(ref[k]-fast['density_A9'][k]))) for k in ('orbital','spin'))
    _,reverse=solve(replace(model,phi=model.phi+np.pi))
    reversal_error=max(float(np.max(np.abs(reverse[p][k]+results[p][k])))
                       for p in results for k in ('orbital','spin'))
    lo=float(g.energies[:,1].max());hi=float(g.energies[:,2].min())
    mus=np.linspace(lo+.2*(hi-lo),hi-.2*(hi-lo),13)
    plateau={p:[integrate(g,float(mu),0,protocol=p)['orbital'][0] for mu in mus] for p in results}
    angles=np.linspace(0,2*np.pi,17)
    angular=[]
    for phi in angles:
        m=replace(model,phi=float(phi));_,rs=solve(m,n=161,mu=.08)
        angular.append({'phi_rad':float(phi),**{
            p:{k:v.imag for k,v in response_si(r,m,.0001,convention='harmonic').items()}
            for p,r in rs.items()}})
    out={
        'schema_version':1,
        'paper':{'doi':'10.1103/xr5v-gkv5','citation':'PRB 114, 245420 (2026)',
                 'basis':'uploaded 17-page publisher PDF; arXiv v1 is not substituted'},
        'acceptance':'MODEL_READY_FOR_EXPLICIT_PROTOCOLS_WITH_OPEN_SOURCE_DISCREPANCY',
        'parameters':{'t0_eV':.08,'t1_eV':1.,'alpha_eV':.8,'exchange_eV':.85,
                      'phi_rad':float(model.phi),'mu_eV':0.,'temperature_K':25.,
                      'drive_energy_eV':.0001,'g':2.,'a_angstrom':3.82,
                      'a_provenance':'explicit illustrative input; in-plane a not specified in uploaded paper',
                      'field_unit_T':model.field_unit_T},
        'convention':'B=B0*cos(wt); chi=j_plus/B0; full phasor=2*chi; omega=(hbar omega)/hbar',
        'benchmark':benchmark,'gap_edges_eV':[lo,hi],
        'convergence':convergence,
        'independent_projector_FD_error':reference_error,
        'neel_reversal_max_error':reversal_error,
        'T0_gap_plateau':{'mu_eV':mus,**plateau},
        'angle_scan':{'mu_eV':.08,'mesh':161,'data':angular},
        'limits':[
            'D10 and A9 orbital results disagree; no automatic correction or merged physical amplitude',
            'CuMnAs doublet-block continuation and fixed orthonormal basis only; not a general Wannier solver',
            'Low frequency, nonresonant intrinsic polarization; no finite-tau transport or Maxwell/circuit model',
            'Zero-temperature metals fail closed; finite-temperature metallic scans need their own mesh convergence',
            'No author raw data, full-figure quantitative certification, or experimental voltage reproduction'
        ]}
    print(json.dumps(out,default=serial,indent=2,allow_nan=False))

if __name__=='__main__':main()
