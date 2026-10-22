"""v65 reproducible finite equilibrium orbital E/B Maxwell experiment, NOT transport."""
import argparse,json
from pathlib import Path
from dataclasses import replace
import numpy as np
from qgmap.equilibrium_orbital_eb_v65 import (
    Spec,benchmark_graph,full_projected_polarization,mixed_static_EB,
    canonical_mixed_static_EB,finite_difference_mixed,finite_difference_canonical)


def generate():
    h,xy,d=benchmark_graph();s=Spec()
    p=full_projected_polarization(xy,s.charge_over_hbar,'y',d)
    k=np.array([[.3,.21+.1j,0,0],[.21-.1j,-.4,.11,0],
                [0,.11,.2,.17+.08j],[0,0,.17-.08j,-.1]],complex)
    data={'version':'v65',
          'scientific_status':'VERIFIED finite equilibrium Maxwell and supplied position-contact algebra; original driven orbital EB UNRESOLVED / source PARTIAL',
          'units':'dimensionless finite 4-orbital model; hbar=1; beta temperature energy, signed q/hbar explicitly set',
          'protocol':'STATIC equilibrium Grand Omega and canonical Helmholtz at finite T, not dynamical Hall/DC transport',
          'parameters':{'mu':s.mu,'temperature':s.temperature,'B':s.field,'charge_over_hbar':s.charge_over_hbar},
          'frames':{},'finite_difference':{},'origin_control':{},'time_reversal_control':{}}
    for frame,K in [('symmetric',None),('landau_y',None),('orbital',k)]:
        data['frames'][frame]={'fixed_mu':mixed_static_EB(h,xy,p,s,frame,K),
                              'fixed_N':canonical_mixed_static_EB(h,xy,p,s,frame,K)}
        data['finite_difference'][frame]={r:float(finite_difference_mixed(h,xy,p,s,frame,K,route=r))
                                         for r in ('magnetic_of_electric','electric_of_magnetic','free_energy')}
        data['finite_difference'][frame].update({'canonical_'+r:float(finite_difference_canonical(h,xy,p,s,frame,K,route=r))
                                         for r in ('magnetic_of_electric','electric_of_magnetic','free_energy')})
    for shift in (0,2):
        data['origin_control'][str(shift)]={
            'fixed_mu':mixed_static_EB(h,xy,p+shift*np.eye(4),s)['mixed_d2F_dEdB'],
            'fixed_N':canonical_mixed_static_EB(h,xy,p+shift*np.eye(4),s)['mixed_d2F_dEdB_fixed_N']}
    tr_h,tr_xy,tr_d=benchmark_graph(False)
    tr_p=full_projected_polarization(tr_xy,s.charge_over_hbar,'y',tr_d)
    for b in (-s.field,0.,s.field):
        data['time_reversal_control'][str(b)]=mixed_static_EB(tr_h,tr_xy,tr_p,replace(s,field=b))['mixed_d2F_dEdB']
    return data


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--output',default='results/frozen/v65_static_equilibrium_orbital_eb.json')
    args=cli.parse_args();path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(generate(),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v65 frozen model evidence',path)
