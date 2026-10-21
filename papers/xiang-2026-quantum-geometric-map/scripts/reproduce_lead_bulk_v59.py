#!/usr/bin/env python3
"""Reproduce v59 1D dispersive-lead Hall device and separate periodic bulk Chern oracle."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.lead_selfenergy_v59 import (
 ChainSpec,reference_device,integrate_leads,finite_B_fd,selfenergy_and_B)
from qgmap.periodic_kubo_v59 import periodic_chern

def evaluate():
    out={'version':'v59',
         'status':'VERIFIED_DEFINED_LEAD_AND_QWZ_ALGORITHMS_CONDITIONAL_DEVICE_NOT_ORBITAL_EB',
         'units':'finite-device qB area/hbar; conductance q²/h; bulk Chern e²/h',
         'lead_reference':[],'periodic_bulk_charge_chern':[],
         'scope':['Not a full orbital magnetoelectric EB response',
                  '1D chain leads are a defined dispersive reservoir, not generic ab-initio leads',
                  'Linear onsite_B_slope is prescribed magnetogating and can violate spinless Onsager',
                  'Periodic QWZ is intrinsic static anomalous charge Hall, separate bulk oracle',
                  'Device Hall voltage and periodic bulk Chern conductance are different observables']}
    spec=ChainSpec()
    for t in (2.,5.,20.,100.):
        H,p,leads=reference_device(lead_hopping=t)
        a=integrate_leads(H,p,leads,spec)
        out['lead_reference'].append({
          'lead_t':t,'coupling_each':leads[0].coupling,
          'gamma_center':selfenergy_and_B(0,0,leads[0])[2],
          'Hall_voltage_per_LR_bias':a['observables']['hall_ratio'],
          'dHall_dB':a['observables']['dHall_dB'],
          'longitudinal_conductance_q2_over_h':a['observables']['longitudinal_conductance_q2_over_h'],
          'quadrature_error':a['quadrature_error'],
          'quadrature_evaluations':a['quadrature_evaluations']})
    H,p,leads=reference_device()
    a=integrate_leads(H,p,leads,spec)
    other=integrate_leads(H,p,leads,spec,gauge='landau_y')
    fd=finite_B_fd(H,p,leads,spec)
    out['independent_benchmark']={
      'hall':a['observables']['hall_ratio'],
      'slope_analytic':a['observables']['dHall_dB'],
      'slope_fd':fd['dHall_dB'],
      'slope_analytic_minus_fd':a['observables']['dHall_dB']-fd['dHall_dB'],
      'max_gauge_T_difference':float(np.max(abs(a['T']-other['T']))),
      'max_gauge_dT_difference':float(np.max(abs(a['dT_dB']-other['dT_dB']))),
      'optical_current_residual':a['current_conservation_residual']}
    H,p,shifted=reference_device(onsite_B_slope=.12)
    a=integrate_leads(H,p,shifted)
    fd=finite_B_fd(H,p,shifted)
    out['field_dependent_selfenergy_test']={
       'onsite_B_slope':.12,'dHall_analytic':a['observables']['dHall_dB'],
       'dHall_five_point':fd['dHall_dB'],
       'note':'Prescribed scalar magnetogating, NOT a microscopic spinless orbital lead'}
    for mass in (-3.,-1.,1.,3.):
        out['periodic_bulk_charge_chern'].append(periodic_chern(mass,31))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    text=json.dumps(evaluate(),indent=2,allow_nan=False)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(text,encoding='utf8')
    else:print(text,end='')
if __name__=='__main__':main()
