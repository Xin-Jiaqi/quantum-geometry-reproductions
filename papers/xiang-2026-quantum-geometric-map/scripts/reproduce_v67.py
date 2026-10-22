"""Independent reproducible noninteracting conserving contact Ward and gate-B benchmark."""
import argparse,json
from dataclasses import replace
from pathlib import Path
from qgmap.reservoir_ward_v67 import (ReservoirSpec,chain_surface,reference,
    scattering,conductance,conductance_B,conductance_E,conductance_EB,
    finite_difference,finite_difference_mixed,local_polarization_bias,magnetic_ward)


def create():
    s=ReservoirSpec()
    out={'version':'v67',
         'status':'VERIFIED noninteracting elastic two-reservoir Ward identities; CONDITIONAL finite open Peierls gate transport; Xiang original EB PARTIAL',
         'observable':'T(mu,E_gate,B), differential linear-reservoir-bias two-terminal conductance; NOT v66 bulk sigma_xy or homogeneous E-driven orbital EB',
         'limits':'one energy mu inside both single-channel lead bands, zero-T differential voltage, 4-site finite device; fixed coherent E/B gates; no disorder/inelastic collisions',
         'parameters':{k:float(getattr(s,k)) for k in s.__dataclass_fields__},
         'surface_g_real':float(chain_surface(s).real),'surface_g_imag':float(chain_surface(s).imag),
         'frames':{},'negative_controls':{},'energy_scan':{}}
    for frame in ('symmetric','landau_y','orbital'):
        op=reference(frame,spec=s);ward=scattering(op);inner=local_polarization_bias(op)
        out['frames'][frame]={
            'conductance_q2_over_h':conductance(op),'dB_conductance':conductance_B(op),
            'dE_gate_conductance':conductance_E(op),'dB_dE_gate_conductance':conductance_EB(op),
            'fd5_dB_conductance':finite_difference(s,frame=frame),
            'fd5_dE_gate_conductance':finite_difference(s,frame=frame,parameter='E'),
            'fd16_dB_dE_gate':finite_difference_mixed(s,frame=frame),
            'scattering_unitarity_residual':ward['unitarity_error'],
            'S_L_to_R_probability':ward['cross_probability'],
            'magnetic_scattering_Ward':magnetic_ward(op),
            'optical_spectral_identity_residual':op['spectral_identity_error'],
            'uncontacted_orbital_continuity_residual':inner['uncontacted_site3_continuity'],
            'uncontacted_orbital_continuity_B_residual':inner['uncontacted_site3_continuity_B'],
            'internal_Pa_current_under_reservoir_bias':inner['internal_pa_bias_current'],
            'dB_internal_Pa_current_under_bias':inner['internal_pa_bias_current_B']}
        if frame!='symmetric':
            out['negative_controls'][frame]={
                'T_if_device_rotated_but_lead_couplings_held_fixed':conductance(reference(frame,transform_contacts=False)),
                'dB_T_if_contact_sigma_B_omitted':conductance_B(op,omit_contact_B=True),
                'dB_dE_T_if_contact_sigma_B_omitted':conductance_EB(op,omit_magnetic_contacts=True),
                'dB_dE_T_if_dB_Pdrive_omitted':conductance_EB(op,omit_position_contact=True)}
    for mu in (-.7,0.,.13,.7):
        out['energy_scan'][str(mu)]=conductance(reference(spec=replace(s,mu=mu)))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',default='results/frozen/v67_conserving_reservoir_ward.json')
    args=p.parse_args();path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(create(),indent=2,sort_keys=True)+'\n',encoding='utf8')
    print('v67 frozen conserving lead benchmark',path)
