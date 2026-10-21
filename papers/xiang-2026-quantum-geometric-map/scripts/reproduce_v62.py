"""Regenerate v62 same-QWZ conditional charge transport and periodic Kubo evidence."""
import json,argparse
from pathlib import Path
from qgmap.qwz_spatial_v62 import spatial_charge_response
from qgmap.qwz_charge_kubo_v62 import qwz_charge_kubo

def create():
    out={'version':'v62','scientific_status':'VERIFIED named methods; CONDITIONAL coherent charge contacts; original orbital E*B UNRESOLVED',
         'temperature_scope':'device 0 K, bulk Kubo may have finite T; never equate unlike protocol',
         'units':'Rxy h/q2; I under unit voltage in q2/h; local J normalized as q2/h',
         'devices':{},'bulk':{}}
    for L,w,m,mu in ((5,5,-1,0),(5,3,-1,0),(7,5,-1,0),(7,3,-1,0),(7,3,1,0),(7,3,3,0),(5,3,-1,1.5)):
        d=spatial_charge_response(L,m,w,mu=mu,eta=1e-4)
        key=f'L{L}_w{w}_m{m:+g}_mu{mu:g}'
        out['devices'][key]={'Rxy_h_over_q2':d['Rxy_h_over_q2'],
           'IL_q2_over_h':d['drive_current_q2_over_h'],
           'x_cut_currents':d['x_cut_currents'].tolist(),
           'max_interior_divergence':float(abs(d['bond_divergence'][1:-1,1:-1]).max()),
           'spectral_sum_error':d['spectral_identity_error'],
           'ldos_edge_00':float(d['spectral_ldos'][0,0]),
           'ldos_center':float(d['spectral_ldos'][L//2,L//2])}
    for m,mu,t in ((-1,0,0),(1,0,0),(3,0,0),(-1,0,.2),(-1,0,1.0),(-1,1.5,0)):
        d=qwz_charge_kubo(m,mu,t,31)
        out['bulk'][f'm{m:+g}_mu{mu:g}_T{t:g}']={'sigma_pauli_q2_over_h':d['sigma_pauli_q2_over_h'],
            'sigma_interband_q2_over_h':d['sigma_interband_q2_over_h'],'min_direct_gap':d['min_direct_gap']}
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/frozen/v62_qwz_bulk_spatial_charge.json')
    a=p.parse_args();f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(create(),sort_keys=True,indent=2)+'\n')
    print(f'v62 frozen evidence -> {f}')
