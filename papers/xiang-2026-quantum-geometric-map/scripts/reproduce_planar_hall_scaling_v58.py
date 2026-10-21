#!/usr/bin/env python3
"""v58 fixed-site-density planar 4-probe coherent Hall size/contact scan.

Adaptive energy integration by default; NO bulk conductivity inference.
Run from ZIP root with PYTHONPATH=src or within PRB paper module.
"""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.planar_hall_scaling_v58 import (
    PlanarSpec,planar_bar,planar_response,adaptive_planar_response)
from qgmap.lead_transport_v57 import LeadSpec


def evaluate(configs,lead,*,integration='adaptive',epsrel=3e-8):
    if integration not in ('adaptive','fixed_gauss_uncertified'):
        raise ValueError('invalid integration protocol')
    rows=[]
    for plan in configs:
        fn=adaptive_planar_response if integration=='adaptive' else planar_response
        x=fn(plan,lead,epsrel=epsrel) if integration=='adaptive' else fn(plan,lead)
        o=x['observables']
        rows.append({'nx':plan.nx,'ny':plan.ny,'sites':x['sample_sites'],
            'aspect_ratio':x['aspect_ratio'],'lattice_site_density':x['lattice_density'],
            'onsite':plan.onsite,'hopping':plan.hopping,'spacing':plan.lattice_spacing,
            'flux_per_plaquette':x['plaquette_flux'],'chemical_potential':lead.mu,
            'coverage_requested':plan.coverage,'coverage_actual':x['contact_site_counts'],
            'gamma_mode':plan.gamma_mode,'gamma_parameter':plan.gamma_per_site,
            'gamma_reservoir_total':x['gamma_total_by_reservoir'],
            'hall_voltage_per_LR_bias':o['hall_ratio'],
            'magnetic_derivative_hall_voltage_per_bias':o['dHall_dB'],
            'hall_resistance_h_over_q2':o['hall_resistance_h_over_q2'],
            'magnetic_derivative_hall_resistance_h_over_q2':o['dHall_resistance_dB'],
            'LR_conductance_q2_over_h':o['longitudinal_conductance_q2_over_h'],
            'longitudinal_resistance_h_over_q2':o.get('longitudinal_resistance_h_over_q2',
                float(1/o['longitudinal_conductance_q2_over_h'])),
            'max_floating_probe_current':o['max_floating_probe_current'],
            'quadrature_error_estimate':x.get('quadrature_error_estimate'),
            'quadrature_integrand_evaluations':x.get('integrand_evaluations'),
            'quadrature_converged':x.get('quadrature_success',False)})
    return {'version':'v58','status':'VERIFIED_MODEL_ALGORITHM_CONDITIONAL_FINITE_COHERENT_DEVICE_NOT_BULK',
        'numerical_protocol':integration,'energy_integration_relative_tolerance':epsrel if integration=='adaptive' else None,
        'lead_spec':lead.__dict__,'results':rows,
        'bulk_limit_verdict':'UNRESOLVED_NO_CONTROLLED_THERMODYNAMIC_OR_INCOHERENT_LIMIT',
        'caveats':['Fixed density means fixed spatial lattice-site density a^-2, NOT fixed carrier density',
                   'Chemical potential and plaquette magnetic flux held fixed across sizes',
                   'Coherent phase-sensitive mesoscopic Hall voltage can be contact/size dependent',
                   'Contact protocol per_site scales total Gamma with contact width; fixed_reservoir_sum does not',
                   'Four-terminal Hall resistance is a device ratio, not an intrinsic bulk tensor',
                   'Energy integration MUST be converged; fixed Gauss mode is uncertified by default',
                   'Spinless square-lattice conventional magnetic Hall response, not original-paper quantum geometry MNHE']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sizes',default='5x5,7x5,7x7,9x5,9x9')
    p.add_argument('--coverage',type=float,default=1.)
    p.add_argument('--gamma-mode',choices=['per_site','fixed_reservoir_sum'],default='per_site')
    p.add_argument('--gamma',type=float,default=.6)
    p.add_argument('--hopping',type=float,default=1.)
    p.add_argument('--a',type=float,default=1.)
    p.add_argument('--onsite',type=float,default=0.)
    p.add_argument('--B',type=float,default=.12)
    p.add_argument('--q-over-hbar',type=float,default=-.87)
    p.add_argument('--mu',type=float,default=.65)
    p.add_argument('--T',type=float,default=.12)
    p.add_argument('--nE',type=int,default=80,help='Only fixed_gauss_uncertified mode; adaptive ignores it')
    p.add_argument('--epsrel',type=float,default=3e-8)
    p.add_argument('--integration',choices=['adaptive','fixed_gauss_uncertified'],default='adaptive')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    dims=[]
    for s in a.sizes.split(','):
        pair=s.strip().lower().split('x')
        if len(pair)!=2:p.error('sizes must be 5x5,7x5,...')
        dims.append(tuple(map(int,pair)))
    specs=[PlanarSpec(nx=x,ny=y,lattice_spacing=a.a,hopping=a.hopping,
        onsite=a.onsite,coverage=a.coverage,gamma_per_site=a.gamma,gamma_mode=a.gamma_mode)
            for x,y in dims]
    lead=LeadSpec(field=a.B,charge_over_hbar=a.q_over_hbar,mu=a.mu,
                  temperature=a.T,broadening=a.gamma,energy_points=a.nE)
    report=evaluate(specs,lead,integration=a.integration,epsrel=a.epsrel)
    data=json.dumps(report,allow_nan=False,indent=2)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(data,encoding='utf8')
    else:print(data,end='')
if __name__=='__main__':main()
