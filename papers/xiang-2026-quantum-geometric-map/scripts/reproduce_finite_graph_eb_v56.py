#!/usr/bin/env python3
"""Independently benchmark finite open-graph AC electric response and d_B sigma.

Units: hbar=1, signed q=q_over_hbar, energy and time reciprocal;
no physical-bulk orbital magnetotransport or SI prediction.
"""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.finite_graph_eb_v56 import Spec,state,mixed_eb,finite_field_derivative,independent_sigma_spectral,equilibrium_current


def encode(z):
    if isinstance(z,complex):return {'real':float(z.real),'imag':float(z.imag)}
    if isinstance(z,(np.float64,np.int64)):return z.item()
    raise TypeError(type(z).__name__)


def evaluate(graph=None,spec=Spec()):
    if graph is None:
        hopping=np.array([[.2,1.,.3,.8],[1.,-.1,1.2,.1],[.3,1.2,.05,.9],[.8,.1,.9,.3]],complex)
        positions=np.array([[0,0],[1.25,0],[1.25,.8],[0,.8]],float)
    else:
        with np.load(graph,allow_pickle=False) as f:
            if not {'hoppings','positions'}<=set(f.files):
                raise ValueError('Graph NPZ must contain hoppings and positions')
            hopping=np.asarray(f['hoppings']);positions=np.asarray(f['positions'])
    data={}
    for gauge in ('symmetric','landau_y'):
        result=mixed_eb(hopping,positions,spec,gauge)
        finite=finite_field_derivative(hopping,positions,spec,gauge,step=1e-4,points=5)
        z=state(hopping,positions,spec,gauge)
        spectral=independent_sigma_spectral(z['H'],z['P_current'],z['P_drive'],spec.mu,spec.temperature,spec.gamma,spec.omega)
        data[gauge]={**result,'finite_B_difference':finite,
                      'analytic_minus_FD':result['mixed_dsigma_dB']-finite,
                      'independent_spectral_sigma':spectral,
                      'equilibrium_current_E0':equilibrium_current(hopping,positions,spec,gauge,0),
                      'equilibrium_current_E0p2':equilibrium_current(hopping,positions,spec,gauge,.2)}
    return {'protocol':'v56_finite_open_graph_E_ac_and_mixed_EB_independent_response',
        'acceptance':'VERIFIED_OPERATOR_AND_DIFFERENTIATION_CONDITIONAL_SCALAR_GAMMA_RESPONSE',
        'scenario':{'mu':spec.mu,'temperature':spec.temperature,'gamma':spec.gamma,
                    'omega':spec.omega,'field_B':spec.field,'charge_over_hbar':spec.q_over_hbar,
                    'electric_driver':'H_E=-E*P_b','probe_axis':spec.probe_axis,'drive_axis':spec.drive_axis,
                    'boundary':'finite open graph','units':'hbar=1, signed q=q_over_hbar'},
        'gauges':data,
        'gauge_mixed_residual':data['symmetric']['mixed_dsigma_dB']-data['landau_y']['mixed_dsigma_dB'],
        'gauge_sigma_residual':data['symmetric']['response_sigma']-data['landau_y']['response_sigma'],
        'limitations':['SCALAR relaxation gamma around magnetic-field-equilibrated f(H_B) is a declared phenomenology',
                       'The finite open-system AC derivative is NOT the intrinsic bulk orbital MNHE/PHE',
                       'Terms individually depend on the B-dependent magnetic gauge',
                       'dc gamma->0 equilibrium current vanishes; finite gamma DC depends on bath model',
                       'Current contacts, equilibrium B derivative and Liouvillian vertex are retained',
                       'No realistic disorder, leads, Landauer current, SI conductivity or thermodynamic limit']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--graph-npz',type=Path,help='optional NPZ; omit to run the fixed four-site reference graph')
    p.add_argument('--mu',type=float,default=.1)
    p.add_argument('--T',type=float,default=.4)
    p.add_argument('--gamma',type=float,default=.08)
    p.add_argument('--omega',type=float,default=.35)
    p.add_argument('--B',type=float,default=.27)
    p.add_argument('--q-over-hbar',type=float,default=-.87)
    p.add_argument('--current-axis',choices=['x','y'],default='x')
    p.add_argument('--electric-axis',choices=['x','y'],default='y')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    result=evaluate(a.graph_npz,Spec(a.mu,a.T,a.gamma,a.omega,a.q_over_hbar,a.B,a.current_axis,a.electric_axis))
    out=json.dumps(result,default=encode,allow_nan=False,indent=2)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(out,encoding='utf-8')
    else:print(out,end='')
if __name__=='__main__':main()
