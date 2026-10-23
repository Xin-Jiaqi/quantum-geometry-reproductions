"""Recompute every v76 physics acceptance number from live operators."""
from dataclasses import asdict
from pathlib import Path
import argparse,json
from qgmap.conserving_two_band_collision_v76 import (
 Spec,CoherentSpec,two_band_populations,population_weighted_spectral,
 bx_fixed_radius_nonelastic,coherent_two_color,coherent_frame_covariance,
 decoherence_nonuniqueness,equilibrium_B_collision_contact,
 mass_zeeman_collision_vertex)

def jsonnum(x):
    if isinstance(x,complex):return {'real':float(x.real),'imag':float(x.imag)}
    if isinstance(x,(float,int,str,bool)):return x
    raise TypeError(str(type(x)))

def generate():
    s=Spec();c=CoherentSpec()
    pop=two_band_populations(s)
    raw=coherent_two_color(c);ind=coherent_two_color(c,'dense')
    gauge=coherent_frame_covariance(c)
    uniqueness=decoherence_nonuniqueness(c)
    bx=equilibrium_B_collision_contact(c,B=1e-4)
    bz=mass_zeeman_collision_vertex(c,delta=1e-5)
    return dict(version='v76',source_science='PARTIAL',
        scope='band-resolved elastic Boltzmann, untilted-coherent CP Lindblad; Bz controlled collision vertex NOT Xiang spin-Bx PHE',
        population_spec=asdict(s),coherent_spec=asdict(c),
        population_total=jsonnum(pop['current']),
        population_per_band={str(b):jsonnum(v['total_current']) for b,v in pop['bands'].items()},
        population_independent_spectral={str(b):jsonnum(population_weighted_spectral(s,b)) for b in (-1,1)},
        fixed_lab_ring_nonelastic=bx_fixed_radius_nonelastic(s),
        coherent={k:jsonnum(raw[k]) for k in ('b_then_e','e_then_b','total','local_harmonic_1','local_mixed_norm','trace_second')},
        independent_dense_local_harmonic=jsonnum(ind['local_harmonic_1']),
        gauge_local_harmonic=jsonnum(gauge['local_harmonic_1']),
        pure_dephasing_no_go={k:jsonnum(uniqueness[k]) for k in
            ('population_L_error','local_harmonic_zero_phi','local_harmonic_phi',
             'local_harmonic_difference','spin_zero_phi','spin_phi','spin_difference')},
        bx_frozen_collision_no_go={k:jsonnum(bx[k]) for k in (
            'analytic_coherent_identity_residual','frozen_B_stationarity_defect',
            'required_collision_contact_norm','required_contact_total_trace','finite_B_derivative_residual')},
        Bz_consistent_mixed={k:jsonnum(bz[k]) for k in (
            'linear_at_B0','analytic_magnetic_derivative','independent_five_point',
            'full_difference_error','frozen_collision_derivative','missing_collision_vertex_error',
            'collision_vertex_drive_norm','equilibrium_collision_ward_residual',
            'field_equilibrium_collision_defect_at_zero')},
        scientific_status='VERIFIED_NAMED_TWO_BAND_MODELS_NOT_ORIGINAL_SOURCE_COMPLETE')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/frozen/v76_conserving_two_band_collision.json')
    args=ap.parse_args();path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(generate(),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('v76 regenerated',path)
