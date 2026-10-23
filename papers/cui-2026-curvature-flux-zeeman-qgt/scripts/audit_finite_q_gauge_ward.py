"""Round27 reproducible finite-q longitudinal gauge identity audit."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.finite_q_gauge_ward import (
 soc_torus,site_phase,bond_currents,divergence,pure_gauge_response,
 linear_liouville,continuity_residual)

t=soc_torus();L=site_phase(t,qx=1,qy=1,offset=.21)
a=pure_gauge_response(t,L,z=.42+.09j,mu=.2,temperature=.3)
ref=linear_liouville(t.H,a['rho'],a['K']+a['Vscalar'],a['z'])
operator_res=1j*(t.H[None]@t.projectors-t.projectors@t.H)+divergence(t,bond_currents(t))
contact_bad=continuity_residual(t,a['dr_vector'],a['K'],a['rho'],a['z'])
contact_good=continuity_residual(t,a['dr_vector'],a['K'],a['rho'],a['z'],a['bond_J_contact'])
report={
'schema':'cui-round27-finite-q-u1-gauge-ward-v1',
'paper_status':'PARTIAL',
'model':'finite 3x3 2D spin-half SOC Peierls torus; no collisions/thermodynamic metallic KME',
'qx_qy':[1,1], 'complex_drive_frequency':[.42,.09],
'operator_continuity_residual_max':float(np.max(abs(operator_res))),
'spectral_vs_orbital_liouville_residual_max':float(np.max(abs(ref-a['dr_gauge']))),
'pure_gauge_density_matrix_residual_max':float(np.max(abs(a['dr_gauge']-a['dr_exact']))),
'scalar_density_abs_max':float(np.max(abs(a['density_scalar']))),
'vector_density_abs_max':float(np.max(abs(a['density_vector']))),
'correct_total_density_abs_max':float(np.max(abs(a['density_gauge']))),
'correct_total_current_abs_max':float(np.max(abs(a['current_gauge']))),
'wrong_missing_contact_current_abs_max':float(np.max(abs(a['current_state_vector']+a['current_state_scalar']))),
'correct_vector_continuity_residual_max':float(np.max(abs(contact_good))),
'wrong_missing_contact_continuity_residual_max':float(np.max(abs(contact_bad))),
'physical_finite_frequency_material_KME_verified':False,
'conserving_scattering_verified':False,
}
for key in ('operator_continuity_residual_max','spectral_vs_orbital_liouville_residual_max',
            'pure_gauge_density_matrix_residual_max','correct_total_density_abs_max',
            'correct_total_current_abs_max','correct_vector_continuity_residual_max'):
 assert report[key] < 1e-10,(key,report[key])
for key in ('wrong_missing_contact_current_abs_max','wrong_missing_contact_continuity_residual_max'):
 assert report[key]>.1,(key,report[key])
frozen=json.loads((ROOT/'results/frozen/round27_finite_q_gauge_ward.json').read_text())
assert report['schema']==frozen['schema'] and frozen['paper_status']=='PARTIAL'
for key,v in report.items():
 if isinstance(v,float):
  assert abs(v-frozen[key])<2e-10, (key,v,frozen[key])
print(json.dumps(report,indent=2,ensure_ascii=False))
