"""v68 frozen noninteracting lead Keldysh EB Ward and same-device observable."""
import argparse,json
from pathlib import Path
from dataclasses import replace
from qgmap.lead_keldysh_eb_v68 import (Spec,response,keldysh_continuity,
                                      finite_difference_B,static_field_FD,static_mixed_FD)
from qgmap.dynamic_orbital_eb_v66 import benchmark_model,dynamic_mixed_eb

def cp(z):return {'real':float(z.real),'imag':float(z.imag)}

def create():
    s=Spec();d={'version':'v68','protocol':'device-only AC length-gauge perturbation -E Pb, microscopic two semiinfinite elastic leads held in fixed Fermi equilibrium',
            'science':'VERIFIED finite frequency Keldysh Dyson Ward and finite-model differentiation; CONDITIONAL finite device and field localization; UNRESOLVED Xiang full bulk orbital E/B and real Wannier embedding',
            'parameters':{'mu':s.mu,'temperature':s.temperature,'B':s.field,'q_over_hbar':s.charge_over_hbar,'omega':s.omega,'lead_hopping':s.lead_hopping,'nodes_per_interval':s.nodes_per_interval},
            'frames':{},'convergence':{},'negative_controls':{},'static_independent':{}}
    for frame in ('symmetric','landau_y','orbital'):
        a=response(s,frame);w=keldysh_continuity(s,frame)
        d['frames'][frame]={'sigma_device':cp(a['sigma_device']),
             'mixed_dB_sigma_device':cp(a['mixed_dB_sigma_device']),
             'electric_current_contact':cp(a['sigma_electric_contact']),
             'mixed_components':{name:cp(v) for name,v in a['mixed_components'].items()},
             'five_point_independent_dB':cp(finite_difference_B(s,frame,step=.005,points=5)),
             'microscopic_ward_error':w['full_ward_error'],
             'uncontacted_site_ward_abs':float(abs(w['uncontacted_site_ward'])),
             'lead_collision_norm':w['lead_collision_norm'],
             'global_number_injection_error':float(abs(w['number_accumulation']-w['number_lead_injection']))}
    for n in (64,384,512,640):
        r=response(replace(s,nodes_per_interval=n))
        d['convergence'][str(n)]={'sigma':cp(r['sigma_device']),'mixed_dB_sigma':cp(r['mixed_dB_sigma_device'])}
    badlead=response(s,'orbital',omit_lead_magnetic_contacts=True)
    badpos=response(s,'orbital',include_all_magnetic_contacts=False)
    h,xy,pa,pb,k=benchmark_model();v66=dynamic_mixed_eb(h,xy,pa,pb)
    d['negative_controls']={'omit_lead_B':cp(badlead['mixed_dB_sigma_device']),
      'omit_position_B':cp(badpos['mixed_dB_sigma_device']),
      'omit_electric_current_contact':cp(response(s,include_electric_contact=False)['sigma_device']),
      'v66_same_operator_but_scalar_gamma_bath':cp(v66['sigma'])}
    dc=replace(s,omega=0.)
    for frame in ('symmetric','landau_y','orbital'):
        a=response(dc,frame)
        d['static_independent'][frame]={'analytic_dE_current':cp(a['sigma_device']),
            'finite_gate_dE_current':cp(static_field_FD(dc,frame,step=.002)),
            'analytic_dB_dE_current':cp(a['mixed_dB_sigma_device']),
            'finite_gate_mixed_dB_dE_current':cp(static_mixed_FD(dc,frame,step=.002))}
    return d

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='results/frozen/v68_same_device_keldysh_eb.json')
    args=parser.parse_args();p=Path(args.output);p.parent.mkdir(exist_ok=True,parents=True)
    p.write_text(json.dumps(create(),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('Saved v68 micro-lead same-device AC EB',p)
