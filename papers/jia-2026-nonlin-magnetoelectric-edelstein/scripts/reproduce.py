"""Reproducible CLI for Jia et al. PRL 136 226302 (2026) method kernels.

Every JSON result names its response *protocol*. Synthetic tight-binding
calculations are NOT experimental CuMnAs/DFT predictions or digitized PRL data.
"""
from __future__ import annotations
import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
import numpy as np

PROJECT=Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:sys.path.insert(0,str(PROJECT))
from nmee_lab import CleanSpec,CuMnAsBGKSpec,evaluate_clean,evaluate_cumnas_bgk,mesh
from nmee_lab.clean_core.models import kane_mele
from nmee_lab.ptbackend.cumnas import Model,points,energies,extrinsic_mu_scan,intrinsic_mu_scan

CITATION='Jia, Xiang, Qiao, Wang, Physical Review Letters 136, 226302 (2026)'

def evaluate(case,n,mu,theta,temperature,tau,a_nm):
    if n<3:raise ValueError('n>=3')
    if case.startswith('fig2'):
        model=kane_mele()
        spec=CleanSpec(alpha='x',gamma='z',beta=0,mu_eV=mu,temperature_K=temperature,branch='printed')
        out=evaluate_clean(model,mesh(model.reciprocal,n),spec)
        # Eq(9) and printed branch are stored together, never conflated.
        return {'case':case,'model':'Kane–Mele Eq.(6)','mesh':[n,n],
                'parameters':{'mu_eV':mu,'temperature_K':temperature,'branch':'both'},
                'protocols':{'prl_printed_eq1':out['parts']['printed'],
                             'static_eq9_full':out['parts']['full_static'],
                             'static_eq9_minus_printed':out['parts']['diagonal_disputed']},
                'normalization':out['normalization'], 'warnings':out['warnings']}
    model=Model(a_nm=a_nm,temperature_K=temperature,tau_fs=tau)
    q=points(n,offset=(.5,.5))  # k/-k exact mesh pairing for symmetry-safe integration
    if case=='fig1b':
        lower,upper,d=energies(q,np.pi/2,model)
        return {'case':case,'model':'CuMnAs Eq.(5)','mesh':[n,n],
                'protocol':'analytic_two_PTx2_doublets',
                'parameters':{'theta_rad':float(np.pi/2),'t_eV':model.t},
                'e_minus_min_eV':float(lower.min()),'e_plus_max_eV':float(upper.max()),
                'min_direct_gap_eV':float(2*d.min()),
                'warnings':['this is a mesh scan, not the full original high-symmetry-path band plot']}
    if case=='fig1c-static-diagnostic':
        # A and B are opposite from PT symmetry, finite-field static susceptibility
        # with frozen occupations, NOT formally the paper's driven DC definition.
        r=intrinsic_mu_scan(q,theta,np.array([mu]),model)[0]
        return {'case':case,'model':'CuMnAs Eq.(5)','mesh':[n,n],
                'protocol':'CuMnAs_frozen_occupation_finiteB_full_static_diagnostic',
                'parameters':{'theta_rad':theta,'mu_eV':mu,'temperature_K':temperature},
                'A_z':float(r[0]),'B_z':float(r[1]),'total_z':float(r[2]),
                'normalization':'historical CuMnAs 2D model units, 1/(nm T V)',
                'warnings':['NOT an independent implementation of PRL Eq.(1) in a degenerate driven dc limit',
                            'Dirac-nearby points can be nonconvergent in the clean static limit']}
    if case=='fig1d-printed':
        # Equation 2 Fermi-surface projector trace, degenerate-safe.
        # p's a and velocity enforce 2D scaling ~1/a.
        v=float(extrinsic_mu_scan(q,theta,np.array([mu]),model,alpha='x',gamma='y',beta=1)[0])
        return {'case':case,'model':'CuMnAs Eq.(5)','mesh':[n,n],
                'protocol':'PRL_Eq2_Fermi_surface_PT_projector_spin_metric_extrinsic_tau',
                'parameters':{'theta_rad':theta,'mu_eV':mu,'temperature_K':temperature,'tau_fs':tau,'a_nm':a_nm},
                'response':v,'normalization':'2D spin magnetic moment in mu_B/nm² per (V/nm*T): 1/(nm V T), Pauli spin convention',
                'warnings':['Fig.1(d) magnitude requires a declared physical lattice length and checked field/spin convention',
                            'uses constant relaxation-time approximation; not a microscopic collision theory']}
    if case=='fig1-bgk':
        s=CuMnAsBGKSpec(theta_rad=theta,mu_eV=mu)
        r=evaluate_cumnas_bgk(s,n=n,params=model)
        return {'case':case,'model':'CuMnAs Eq.(5)','mesh':[n,n],**r}
    raise ValueError('unsupported case')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=['fig2-printed-vs-static','fig1b','fig1c-static-diagnostic','fig1d-printed','fig1-bgk'],required=True)
    parser.add_argument('--mesh',type=int,default=41)
    parser.add_argument('--mu',type=float,default=.025)
    parser.add_argument('--theta',type=float,default=.15,help='Néel angle radians for Fig.1')
    parser.add_argument('--temperature',type=float,default=10.,help='temperature K; use 100 for Fig.1')
    parser.add_argument('--tau',type=float,default=50.,help='momentum lifetime fs for CuMnAs')
    parser.add_argument('--a-nm',type=float,default=1.,help='CuMnAs lattice length nm, explicit model convention')
    parser.add_argument('--output',type=Path)
    a=parser.parse_args(argv)
    r={'citation':CITATION,'result':evaluate(a.case,a.mesh,a.mu,a.theta,a.temperature,a.tau,a.a_nm)}
    txt=json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(txt+'\n',encoding='utf-8')
    print(txt)

if __name__=='__main__':main()
