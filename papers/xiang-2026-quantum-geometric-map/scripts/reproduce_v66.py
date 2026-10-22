"""Frozen complete seven-contact finite AC orbital-EB scalar-gamma benchmark."""
import json,argparse
from pathlib import Path
from dataclasses import replace
from qgmap.dynamic_orbital_eb_v66 import (
    Spec,benchmark_model,ac_state,ac_spectral_independent,dynamic_mixed_eb,
    full_field_finite_difference,independent_dc_mixed_fd)


def pair(z): return {'real':float(z.real),'imag':float(z.imag)}


def create():
    h,xy,pa,pb,K=benchmark_model();s=Spec()
    out={'version':'v66','model':'explicit finite four-orbital Peierls and supplied noncommuting projected P_x/P_y',
         'scientific_status':'VERIFIED seven-contact finite-model derivative; CONDITIONAL scalar gamma bath; original Xiang orbital EB UNRESOLVED, source PARTIAL',
         'units':'dimensionless hbar=1; complex finite-frequency linear-E response and derivative B',
         'parameters':{'mu':s.mu,'temperature':s.temperature,'field':s.field,'charge_over_hbar':s.charge_over_hbar,'gamma':s.gamma,'omega':s.omega},
         'frames':{},'DC_non_equilibrium_oracle':{},'static_contact_control':{}}
    for frame,k in (('symmetric',None),('landau_y',None),('orbital',K)):
        r=dynamic_mixed_eb(h,xy,pa,pb,s,frame,k)
        fd=full_field_finite_difference(h,xy,pa,pb,s,frame,k)
        spectral=ac_spectral_independent(h,xy,pa,pb,s,frame,k)
        out['frames'][frame]={'sigma':pair(r['sigma']),'mixed_derivative':pair(r['mixed_dsigma_dB']),
             'components':{n:pair(v) for n,v in r['components'].items()},
             'five_point_B_derivative':pair(fd),'spectral_independent_sigma':pair(spectral),
             'analytic_vs_FD_error':float(abs(r['mixed_dsigma_dB']-fd))}
        dc=replace(s,omega=0)
        analytical=dynamic_mixed_eb(h,xy,pa,pb,dc,frame,k)['mixed_dsigma_dB']
        numerical=independent_dc_mixed_fd(h,xy,pa,pb,dc,frame,k)
        out['DC_non_equilibrium_oracle'][frame]={'analytical':pair(analytical),
            'independent_finite_E_B_superoperator':pair(numerical),
            'absolute_error':float(abs(analytical-numerical))}
    near=replace(s,omega=0,gamma=1e-6)
    actual=ac_state(h,xy,pa,pb,near)
    wrong=ac_state(h,xy,pa,pb,near,electric_contact=False)
    out['static_contact_control']={'gamma':near.gamma,'with_E_current_contact':pair(actual['sigma']),
       'without_E_current_contact':pair(wrong['sigma'])}
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results/frozen/v66_dynamic_orbital_mixed_EB.json')
    args=parser.parse_args();f=Path(args.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(create(),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v66 finite dynamic orbital EB frozen benchmark:',f)
