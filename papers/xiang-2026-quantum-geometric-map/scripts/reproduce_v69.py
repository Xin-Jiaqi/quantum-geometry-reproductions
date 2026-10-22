"""Reproduce v69 finite-lead Lehmann benchmarks against the v68 continuum Keldysh result."""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.finite_lead_lehmann_v69 import (Spec,lehmann_response,static_isothermal_response,
    static_finite_difference,lehmann_continuity_ward,finite_difference_B)
from qgmap.lead_keldysh_eb_v68 import Spec as KSpec,response as keldysh_response


def pair(z): return {'real':float(np.real(z)),'imag':float(np.imag(z))}


def create(include_production=False):
    kr=keldysh_response(KSpec(nodes_per_interval=640))
    frames={}
    s=Spec(length=60,eta=.03)
    for frame in ('symmetric','landau_y','orbital'):
        r=lehmann_response(s,frame)
        frames[frame]={'sigma_device':pair(r['sigma_device']),
            'sigma_noncontact':pair(r['sigma_noncontact']),
            'sigma_electric_contact':pair(r['sigma_electric_contact'])}
    medium=lehmann_response(Spec(length=250,eta=.02))['sigma_device']
    long=lehmann_response(Spec(length=350,eta=.02))['sigma_device']
    intermediate=lehmann_response(Spec(length=500,eta=.01))['sigma_device']
    ss=Spec(length=50,eta=.02,omega=0.)
    st=static_isothermal_response(ss)['sigma_device']
    stfd=static_finite_difference(ss,step=1e-4)
    ward=lehmann_continuity_ward(Spec(length=40,eta=.03))
    out={'version':'v69','model':'closed finite Hermitian four-orbital Peierls device plus two finite tight-binding leads; eigenstate-pair Lehmann sum; no self-energy integral',
        'scientific_status':'VERIFIED finite-system Lehmann and isothermal identities; CONDITIONAL finite-length eta regularization; original Xiang bulk orbital EB remains PARTIAL',
        'parameters':{'mu':.13,'temperature':.37,'field':.23,'charge_over_hbar':-.84,'omega':.29,
                      'lead_hopping':1.1,'lead_onsite':0.},
        'keldysh_reference':{'sigma_device':pair(kr['sigma_device']),
            'mixed_dB_sigma_device':pair(kr['mixed_dB_sigma_device']),
            'nodes_per_interval':640},
        'quick':{'frames_length60_eta003':frames,
            'length_broadening_sequence':[
                {'length':250,'eta':.02,'sigma_device':pair(medium),'absolute_error_to_keldysh':float(abs(medium-kr['sigma_device']))},
                {'length':350,'eta':.02,'sigma_device':pair(long),'absolute_error_to_keldysh':float(abs(long-kr['sigma_device']))},
                {'length':500,'eta':.01,'sigma_device':pair(intermediate),'absolute_error_to_keldysh':float(abs(intermediate-kr['sigma_device']))}],
            'static_isothermal_sigma':pair(st),'static_full_E_finite_difference':pair(stfd),
            'static_abs_error':float(abs(st-stfd)),
            'closed_system_ward':{'ward_error':ward['ward_error'],
                'total_charge_response':pair(ward['total_charge_response']),
                'device_charge_response':pair(ward['device_charge_response']),
                'lead_charge_response':pair(ward['lead_charge_response']),
                'charge_partition_error':ward['charge_partition_error']}},
        'limit_protocol':'Jointly increase finite lead length L and decrease retarded response broadening eta; eta is a spectral resolution, not a phenomenological reservoir bath.'}
    if include_production:
        sprod=Spec(length=1000,eta=.005)
        prod=lehmann_response(sprod)['sigma_device']
        out['production']={'length':1000,'eta':.005,'sigma_device':pair(prod),
            'absolute_error_to_keldysh':float(abs(prod-kr['sigma_device'])),
            'interpretation':'Long finite-lead Lehmann benchmark; residual difference is finite-L/eta error and is not declared zero.'}
        out['production']['five_point_B_derivative_length350_eta002']=pair(
            finite_difference_B(Spec(length=350,eta=.02),step=.003,points=5))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='results/frozen/v69_finite_lead_lehmann.json')
    p.add_argument('--quick',action='store_true',help='omit L=1000 production run')
    a=p.parse_args();path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(create(not a.quick),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v69 finite-lead Lehmann benchmark written:',path)
