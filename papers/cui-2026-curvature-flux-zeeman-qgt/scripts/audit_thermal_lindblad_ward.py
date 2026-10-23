"""Round29 bounded thermal-vs-local CP bath audit and Ward obstruction.

This audit first records measured values; a later CI-checked frozen JSON
is mandatory for acceptance, not inferred from the script's printout.
"""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.thermal_lindblad_ward import (
 two_site_soc,canonical_gibbs,local_lindblad,thermal_transition_jumps,
 dissipator,liouvillian,continuity_terms,bath_gauge_vertex,
 gauge_rotate,evolved_rk4,evolved_dense,bond_current)
H,P,A=two_site_soc();beta=.8;gamma=.16
rho=canonical_gibbs(H,beta)
local=local_lindblad(P,gamma=.18)
bath=thermal_transition_jumps(H,A,beta=beta,gamma=gamma)
jumps=bath['jumps'];rates=bath['rates'];e=bath['energies']
p=np.exp(-beta*(e-e.min()));p/=p.sum()
X=P[0]/2
G=np.diag([.47,.47,-.33,-.33]).astype(complex)
phi=.39
Uphase=lambda Y:gauge_rotate(Y,G,phi)
rot_jumps=np.array([Uphase(j) for j in jumps])
dr=1j*(G@rho-rho@G)
Psource=np.array([np.trace(Pi@dissipator(X,jumps)) for Pi in P])
localterms=continuity_terms(H,X,P,local)
thermalterms=continuity_terms(H,X,P,jumps)
joint_residual=dissipator(Uphase(X),rot_jumps)-Uphase(dissipator(X,jumps))
mismatch=dissipator(Uphase(X),jumps)-Uphase(dissipator(X,jumps))
eps=2e-5
def driven(phi):
 J=np.array([gauge_rotate(j,G,phi) for j in jumps])
 return dissipator(rho,J)
fd=(driven(eps)-driven(-eps))/(2*eps)
vertex=bath_gauge_vertex(rho,jumps,G)
rk=evolved_rk4(H,X,jumps,dt=.01,steps=150)
exact=evolved_dense(H,X,jumps,time=1.5)
report={
 'schema':'cui-round29-lindblad-local-charge-vs-thermal-balance-v1',
 'paper_status':'PARTIAL',
 'model':'synthetic two-site spin-half SOC single-particle canonical Gibbs, not metal KME',
 'beta':beta,'gamma_thermal':gamma,'gamma_local':.18,
 'thermal_pair_detailed_balance_max':float(np.max(abs(rates*p[None,:]-rates.T*p[:,None]))),
 'thermal_gibbs_stationarity_abs_max':float(np.max(abs(dissipator(rho,jumps)))),
 'local_dephasing_gibbs_stationarity_violation_max':float(np.max(abs(dissipator(rho,local)))),
 'local_dephasing_gibbs_energy_drift':float(np.trace(H@dissipator(rho,local)).real),
 'thermal_bath_local_charge_source_abs_max':float(np.max(abs(Psource))),
 'thermal_bath_net_particle_source_abs':float(abs(Psource.sum())),
 'local_bath_charge_continuity_residual_max':float(np.max(abs(localterms['local_residual']))),
 'thermal_bath_charge_continuity_with_source_residual_max':float(np.max(abs(thermalterms['local_residual']))),
 'thermal_bath_wrong_omitted_source_residual_max':float(np.max(abs(thermalterms['density_rate']+thermalterms['divergence']))),
 'thermal_joint_local_U1_covariance_residual_max':float(np.max(abs(joint_residual))),
 'wrong_stale_jump_gauge_error_max':float(np.max(abs(mismatch))),
 'bath_vertex_finite_difference_residual_max':float(np.max(abs(fd-vertex))),
 'wrong_omitted_bath_vertex_gauge_error_max':float(np.max(abs(dissipator(dr,jumps)))),
 'complete_bath_pure_gauge_cancel_residual_max':float(np.max(abs(dissipator(dr,jumps)+vertex))),
 'exact_dense_vs_RK4_max':float(np.max(abs(exact-rk))),
 'rk4_min_eigenvalue':float(np.min(np.linalg.eigvalsh((rk+rk.conj().T)/2))),
 'microscopic_scattering_derived':False,
 'local_electron_charge_ward_for_thermal_bath_without_source':False,
 'full_material_KME':False
}
for key in ('thermal_pair_detailed_balance_max','thermal_gibbs_stationarity_abs_max',
 'thermal_bath_net_particle_source_abs','local_bath_charge_continuity_residual_max',
 'thermal_bath_charge_continuity_with_source_residual_max',
 'thermal_joint_local_U1_covariance_residual_max',
 'bath_vertex_finite_difference_residual_max',
 'complete_bath_pure_gauge_cancel_residual_max','exact_dense_vs_RK4_max'):
 assert report[key]<2e-8,(key,report[key])
for key in ('local_dephasing_gibbs_stationarity_violation_max',
 'thermal_bath_local_charge_source_abs_max',
 'thermal_bath_wrong_omitted_source_residual_max',
 'wrong_stale_jump_gauge_error_max',
 'wrong_omitted_bath_vertex_gauge_error_max'):
 assert report[key]>1e-3,(key,report[key])
assert report['rk4_min_eigenvalue']>-1e-9
ref_file=ROOT/'results/frozen/round29_thermal_lindblad_ward_audit.json'
ref=json.loads(ref_file.read_text())
assert ref['schema']==report['schema'] and ref['paper_status']=='PARTIAL'
for k,v in report.items():
 if isinstance(v,float):
  assert abs(v-ref[k])<3e-10,(k,v,ref[k])
print(json.dumps(report,indent=2))
