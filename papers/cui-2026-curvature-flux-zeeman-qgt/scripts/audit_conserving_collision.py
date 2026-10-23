"""Round28 independent finite-q charge-conserving toy collision audit.

Frozen data becomes an acceptance gate in a subsequent reviewed commit.
"""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.finite_q_gauge_ward import soc_torus,site_phase
from flux_lab.conserving_collision import (
    collision_response,pure_gauge_drive,zeeman_vertex)

t=soc_torus();L=site_phase(t,qx=1,qy=1,offset=.21)
z=.42+.13j;gamma=.12
K,Vs=pure_gauge_drive(t,L,z)
Vz=zeeman_vertex(t,np.diagonal(L).real[::2],spin_axis='y',kappa=.7)
pure=collision_response(t,K+Vs,G=L,z=z,gamma=gamma)
scalar=collision_response(t,.17*L,z=z,gamma=gamma)
spin=collision_response(t,Vz,z=z,gamma=gamma)
other=collision_response(t,Vz+K+Vs,G=L,z=z,gamma=gamma)
wrong_rta=collision_response(t,.17*L,z=z,gamma=gamma,
                             conserve_local_charge=False,gauge_reference=False)
wrong_reference=collision_response(t,K+Vs,G=L,z=z,gamma=gamma,gauge_reference=False)
dense=collision_response(t,Vz+K+Vs,G=L,z=z,gamma=gamma,method='liouville')
dr_expected=1j*(L@spin['rho']-spin['rho']@L)
out={
 'schema':'cui-round28-local-charge-conserving-gauge-covariant-linear-collision-v1',
 'paper_status':'PARTIAL',
 'model':'3x3 SOC single-particle finite torus; affine linearized site dephasing; not a CP thermal bath',
 'gamma':gamma,'drive_z':[z.real,z.imag],
 'pure_gauge_density_abs_max':float(np.max(abs(pure['density']))),
 'pure_gauge_current_abs_max':float(np.max(abs(pure['current']))),
 'pure_gauge_density_matrix_oracle_residual_max':float(np.max(abs(pure['dr']-dr_expected))),
 'scalar_density_abs_max':float(np.max(abs(scalar['density']))),
 'scalar_current_abs_max':float(np.max(abs(scalar['current']))),
 'scalar_charge_ward_residual_max':float(np.max(abs(scalar['ward']))),
 'zeeman_charge_ward_residual_max':float(np.max(abs(spin['ward']))),
 'spin_plus_em_gauge_covariance_density_matrix_error_max':float(np.max(abs(other['dr']-spin['dr']-dr_expected))),
 'spin_plus_em_gauge_covariance_current_error_max':float(np.max(abs(other['current']-spin['current']))),
 'collision_site_charge_residual_max':float(np.max(abs(spin['collision_site']))),
 'dense_vs_iterative_liouville_residual_max':float(np.max(abs(dense['dr']-other['dr']))),
 'naive_rta_charge_ward_violation_max':float(np.max(abs(wrong_rta['ward']))),
 'missing_covariant_collision_reference_pure_gauge_current_max':float(np.max(abs(wrong_reference['current']))),
 'conserving_micro_bath_claim':False,
 'real_material_zeeman_KME_claim':False,
}
# Pivotal audit: separate charge conservation from gauge invariance.
for key in ('pure_gauge_density_abs_max','pure_gauge_current_abs_max',
            'scalar_charge_ward_residual_max','zeeman_charge_ward_residual_max',
            'spin_plus_em_gauge_covariance_density_matrix_error_max',
            'spin_plus_em_gauge_covariance_current_error_max',
            'collision_site_charge_residual_max',
            'dense_vs_iterative_liouville_residual_max'):
 assert out[key]<2e-10,(key,out[key])
for key in ('naive_rta_charge_ward_violation_max',
            'missing_covariant_collision_reference_pure_gauge_current_max'):
 assert out[key]>1e-3,(key,out[key])
# Reference file is added after the first independent remote run.
ref_path=ROOT/'results/frozen/round28_conserving_collision_audit.json'
if ref_path.exists():
    ref=json.loads(ref_path.read_text())
    assert ref['schema']==out['schema'] and ref['paper_status']=='PARTIAL'
    for key,value in out.items():
        if isinstance(value,float):
            assert abs(value-ref[key])<2e-10,(key,value,ref[key])
print(json.dumps(out,indent=2,ensure_ascii=False))
