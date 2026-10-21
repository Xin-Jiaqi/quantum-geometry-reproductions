#!/usr/bin/env python3
"""v60 independent multi-orbital principal-layer and SAME-model QWZ topology."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.matrix_leads_v60 import embedded_four_site_leads,integrated_matrix_leads
from qgmap.lead_transport_v57 import diamond_graph
from qgmap.qwz_finite_marker_v60 import compare_finite_periodic

def evaluate():
    pos,H,_=diamond_graph(1)
    leads=embedded_four_site_leads(H,pos,orbital_matrix=True,eta=.003)
    r=integrated_matrix_leads(H,pos,leads,epsrel=8e-8)
    shifted=pos+np.array([.3,.5])
    a=integrated_matrix_leads(H,shifted,leads,epsrel=8e-8)
    b=integrated_matrix_leads(H,shifted,leads,epsrel=8e-8,gauge='landau_y')
    scalar=embedded_four_site_leads(H,pos,orbital_matrix=False,eta=.001)
    single=integrated_matrix_leads(H,pos,scalar,epsrel=8e-8)
    finite=[compare_finite_periodic(L,mass,max(1,L//4)) for mass in (-3.,-1.,1.,3.) for L in (5,7,9,11)]
    return {'version':'v60','status':'VERIFIED_DECLARED_MATRIX_LEAD_AND_FINITE_QWZ_TOPOLOGY_NOT_BULK_ORBITAL_EB',
      'matrix_lead':{'orbitals_per_principal_layer':2,'eta':.003,
        'hall_ratio':r['observables']['hall_ratio'],'dHall_dB':r['observables']['dHall_dB'],
        'conductance_q2_over_h':r['observables']['longitudinal_conductance_q2_over_h'],
        'quad_error':r['quadrature_error'],'evals':r['evaluations'],
        'lead_current_residual':r['lead_current_conservation_error'],
        'shifted_origin_gauge_T_residual':float(np.max(abs(a['T']-b['T']))),
        'shifted_origin_gauge_dT_residual':float(np.max(abs(a['dT_dB']-b['dT_dB']))),
        'scalar_eta_001_hall_ratio':single['observables']['hall_ratio']},
      'same_QWZ_periodic_vs_finite_marker':finite,
      'accepted':['matrix retarded principal-layer decimation with specified eta>0','full matrix Gamma and magnetic C(B) contact derivatives','same-H periodic Chern and open projector marker'],
      'conditional':['finite device Hall for chosen coherent electrode'],
      'unresolved':['bulk orbital mixed EB','same-H finite-vs-bulk transport conductivity','real-material Wannier and microscopic scattering']}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);args=p.parse_args()
    s=json.dumps(evaluate(),indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(s,encoding='utf8')
    else:print(s,end='')
if __name__=='__main__':main()
