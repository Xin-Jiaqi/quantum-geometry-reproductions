"""Round30 frozen-auditable *many-electron* Fock CAR/KMS/gauge/Zeeman pilot."""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.fermionic_fock_bath import (
 fermion_ops,fock_soc_hubbard,grand_gibbs,thermal_sector_jumps,
 local_charge_balance,unitary_phase,bath_contact,sector_gibbs_weights)
from flux_lab.thermal_lindblad_ward import dissipator,liouvillian,master_rhs

beta=.83;mu=.16;gamma=.14
d=fock_soc_hubbard()
H=d['H'];N=d['N'];A=d['coupling'];G=d['gauge']
rho=grand_gibbs(H,N,beta=beta,mu=mu)
t=thermal_sector_jumps(H,N,A,beta=beta,gamma=gamma)
J=t['jumps'];r=t['rates'];E=t['energies'];Ns=t['numbers']
p=sector_gibbs_weights(E,Ns,beta=beta,mu=mu)
fock_state=np.zeros((16,16),complex);fock_state[0b0011,0b0011]=1.
balance=local_charge_balance(d,fock_state,J)
rng=np.random.default_rng(30)
z=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
X=z@z.conj().T;X/=np.trace(X)
D=dissipator(X,J)
K=(liouvillian(H,J)@X.reshape(-1,order='F')).reshape((16,16),order='F')
Ldirect=master_rhs(H,X,J)
phi=.41
R=lambda Y:unitary_phase(Y,G,phi)
JR=np.array([R(l) for l in J])
gauged_err=np.max(abs(dissipator(R(X),JR)-R(dissipator(X,J))))
stale_err=np.max(abs(dissipator(R(X),J)-R(dissipator(X,J))))
dr=1j*(G@rho-rho@G)
bcontact=bath_contact(rho,J,G)
eps=2e-5
def Dp(sign):
 js=np.array([unitary_phase(l,G,sign*eps) for l in J])
 return dissipator(rho,js)
fdg=(Dp(1)-Dp(-1))/(2*eps)
b_eps=3e-5; V=d['zeeman']
def eval_b(sign):
 hb=H+sign*b_eps*V
 rb=grand_gibbs(hb,N,beta=beta,mu=mu)
 jb=thermal_sector_jumps(hb,N,A,beta=beta,gamma=gamma)['jumps']
 return rb,jb
rp,jp=eval_b(1);rm,jm=eval_b(-1)
drb=(rp-rm)/(2*b_eps)
Bvertex=(dissipator(rho,jp)-dissipator(rho,jm))/(2*b_eps)
Bstate=dissipator(drb,J)
bare=fock_soc_hubbard(interaction=0,cross=0)
rb=grand_gibbs(bare['H'],bare['N'],beta=beta,mu=mu)
occ=np.array([np.trace(o@rho).real for o in d['occupations']])
occ0=np.array([np.trace(o@rb).real for o in bare['occupations']])
c,cd=fermion_ops()
report={
 "schema":"cui-round30-interacting-four-mode-fermionic-kms-ward-v1",
 "paper_status":"PARTIAL",
 "model":"2 sites x true spin-1/2; 4 fermion modes; 16-state interacting Fock, chosen secular KMS jumps",
 "beta":beta,"mu":mu,"gamma":gamma,
 "hubbard_U":d['interaction'],"density_V":d['cross'],
 "fock_dimension":16,"number_sector_sizes":[int(np.sum(Ns==i)) for i in range(5)],
 "number_conserving_kms_jump_count":int(len(J)),
 "car_max_residual":float(max(np.max(abs(c[i]@cd[j]+cd[j]@c[i]-(np.eye(16) if i==j else 0))) for i in range(4) for j in range(4))),
 "pauli_blocked_hop_amplitude":float(np.linalg.norm(cd[1]@c[0]@np.eye(16)[:,3])),
 "fermionic_occupation_min":float(np.min(occ)),
 "fermionic_occupation_max":float(np.max(occ)),
 "interaction_vs_noninteracting_occupation_max":float(np.max(abs(occ-occ0))),
 "grand_KMS_pair_residual_max":float(np.max(abs(r*p[None,:]-r.T*p[:,None]))),
 "grand_gibbs_GKSL_stationarity_max":float(np.max(abs(dissipator(rho,J)))),
 "arbitrary_density_total_charge_bath_source_abs":float(abs(np.trace(N@D))),
 "arbitrary_density_total_charge_variance_source_abs":float(abs(np.trace(N@N@D))),
 "local_bath_source_abs_max":float(np.max(abs(balance['bath_source']))),
 "sitewise_continuity_with_bath_source_max":float(np.max(abs(balance['continuity']))),
 "wrong_omitted_bath_source_local_continuity_max":float(np.max(abs(balance['density_rate']+np.array([balance['current'],-balance['current']])))),
 "orbital_liouvillian_vs_direct_GKSL_max":float(np.max(abs(K-Ldirect))),
 "joint_jump_gauge_covariance_max":float(gauged_err),
 "wrong_stale_jump_gauge_max":float(stale_err),
 "finite_gauge_bath_vertex_difference_max":float(np.max(abs(fdg-bcontact))),
 "missing_gauge_bath_vertex_error_max":float(np.max(abs(dissipator(dr,J)))),
 "complete_gauge_bath_vertex_residual_max":float(np.max(abs(dissipator(dr,J)+bcontact))),
 "magnetic_thermal_bath_vertex_residual_max":float(np.max(abs(Bstate+Bvertex))),
 "missing_magnetic_bath_vertex_error_max":float(np.max(abs(Bstate))),
 "chemical_potential_relaxation_proved":False,
 "microscopic_electron_phonon_vertex_proved":False,
 "source_paper_material_KME_verified":False
}
for key in ("car_max_residual","pauli_blocked_hop_amplitude","grand_KMS_pair_residual_max",
 "grand_gibbs_GKSL_stationarity_max","arbitrary_density_total_charge_bath_source_abs",
 "arbitrary_density_total_charge_variance_source_abs",
 "sitewise_continuity_with_bath_source_max","orbital_liouvillian_vs_direct_GKSL_max",
 "joint_jump_gauge_covariance_max","finite_gauge_bath_vertex_difference_max",
 "complete_gauge_bath_vertex_residual_max","magnetic_thermal_bath_vertex_residual_max"):
 assert report[key]<2e-7,(key,report[key])
for key in ("interaction_vs_noninteracting_occupation_max",
 "local_bath_source_abs_max","wrong_omitted_bath_source_local_continuity_max",
 "wrong_stale_jump_gauge_max","missing_gauge_bath_vertex_error_max",
 "missing_magnetic_bath_vertex_error_max"):
 assert report[key]>1e-5,(key,report[key])
assert report['number_sector_sizes']==[1,4,6,4,1]
assert report['fermionic_occupation_min']>0 and report['fermionic_occupation_max']<1
freeze=ROOT/'results/frozen/round30_interacting_fock_bath_audit.json'
expected=json.loads(freeze.read_text())
assert expected['schema']==report['schema'] and expected['paper_status']=="PARTIAL"
for key,value in report.items():
 if isinstance(value,float):
  assert abs(value-expected[key])<2e-9,(key,value,expected[key])
print(json.dumps(report,indent=2,ensure_ascii=False))
