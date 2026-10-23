"""Frozen charge-conserving Dirac isofrequency scattering and source-sign gates."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.elastic_spin_collision_v75 import Spec,collision_data,eig_physics,steady,sign_invariance

def to_json(a):
    if isinstance(a,(np.floating,np.integer)):return a.item()
    if isinstance(a,(complex,np.complexfloating)):return {'real':float(np.real(a)),'imag':float(np.imag(a))}
    if isinstance(a,np.ndarray):return [to_json(x) for x in a.tolist()]
    if isinstance(a,tuple):return [to_json(x) for x in a]
    if isinstance(a,dict):return {k:to_json(v) for k,v in a.items()}
    return a

def generate():
    specs={'untilted_scalar':Spec(beta=0,spin_overlap=0),
           'untilted_dirac_spinor':Spec(beta=0,spin_overlap=1),
           'tilted_dirac_spinor':Spec(beta=.125,spin_overlap=1),
           'stronger_tilt_dirac_spinor':Spec(beta=.5,spin_overlap=1)}
    out={}
    for key,s in specs.items():
        e=eig_physics(s);sol=steady(s)
        out[key]={'beta':s.beta,'overlap':s.spin_overlap,'gamma':s.gamma,
                   'tau_current':to_json(sol['tau_current']),
                   'gamma_times_tau_current':to_json(s.gamma*sol['tau_current']),
                   'spectral_gap':e['spectral_gap'],
                   'charge_conservation_error':e['number_conservation_error'],
                   'detailed_balance_error':e['detailed_balance_error'],
                   'energy_shell_error':e['energy_shell_error'],
                   'equal_radius_energy_spread':e['equal_radius_energy_spread']}
    rng=np.random.default_rng(734)
    a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));r=(a+a.conj().T)/2
    b=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));sp=(b+b.conj().T)/2
    sign=sign_invariance(np.array([-.8,-.21,.51,1.2]),r,sp,np.array([.3,-.11,.19,-.25]))
    return {'version':'v75','scientific_status':'PARTIAL',
            'scope':'elastic intra-band energy shell transport and literal source sign no-go; not full original spin PHE or orbital MNHE',
            'collision_models':out,
            'sign_audit':{'eq20_direct':to_json(sign['literal']['direct_eq20']),
                          'printed_eq31':to_json(sign['literal']['printed_eq31']),
                          'relative_sum':to_json(sign['relative_sum'])}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/frozen/v75_elastic_collision_and_source_sign.json')
    a=p.parse_args();f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(generate(),sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('generated',f)
