"""Round31 independently audited bosonic oscillator-spectrum/Fock dynamic Ward.

First independent remote run records raw metrics. Later checked snapshot mandatory.
Units are finite synthetic, z=drive omega+i eta, not macroscopic metallic KME.
"""
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from flux_lab.fermionic_fock_bath import (
    fock_soc_hubbard,grand_gibbs,local_charge_balance,unitary_phase)
from flux_lab.thermal_lindblad_ward import dissipator
from flux_lab.bosonic_davies_vertex import (
    finite_oscillator_correlator,bose_davies_fock_jumps,
    ohmic_spectral_density,oscillator_fourier_rate,
    bath_magnetic_derivative,linear_response,bath_contact)

beta=.83;mu=.16;g2=.045;cutoff=2.4
data=fock_soc_hubbard()
H=data['H'];N=data['N'];G=data['gauge']
spec=bose_davies_fock_jumps(H,N,data['coupling'],beta=beta,g2=g2,cutoff=cutoff)
rho=grand_gibbs(H,N,beta=beta,mu=mu)
E=spec['energies'];ns=spec['numbers'];W=spec['rates'];jumps=spec['jumps']
p=np.exp(-beta*(E-mu*ns));p/=p.sum()
state=np.zeros((16,16),complex);state[3,3]=1
local=local_charge_balance(data,state,jumps)
os=finite_oscillator_correlator(1.25,beta=beta,levels=48)
B=os['b']+os['b'].conj().T;U=np.diag(np.exp(1j*.37*np.diag(os['H']).real))
corr_mat=np.trace(os['rho']@(U@B@U.conj().T)@B)
corr_pair=os['emission']*np.exp(-1j*1.25*.37)+os['absorption']*np.exp(1j*1.25*.37)
mag=bath_magnetic_derivative(data,beta=beta,g2=g2,cutoff=cutoff,mu=mu)
g=linear_response(data,mode='gauge',beta=beta,mu=mu,g2=g2,cutoff=cutoff)
gwrong=linear_response(data,mode='gauge',beta=beta,mu=mu,g2=g2,cutoff=cutoff,include_bath_contact=False)
b=linear_response(data,mode='magnetic',beta=beta,mu=mu,g2=g2,cutoff=cutoff)
bwrong=linear_response(data,mode='magnetic',beta=beta,mu=mu,g2=g2,cutoff=cutoff,include_bath_contact=False)
exact_g=1j*(G@rho-rho@G)
eps=2e-5
def dp(s):
    js=np.array([unitary_phase(j,G,s*eps) for j in jumps])
    return dissipator(rho,js)
contact_fd=(dp(1)-dp(-1))/(2*eps)
phi=.37
R=lambda X:unitary_phase(X,G,phi)
Rjs=np.array([R(j) for j in jumps])
X=.7*state+.3*rho
report=dict(
 schema='cui-round31-bosonic-correlation-davies-finite-z-fock-ward-v1',
 paper_status='PARTIAL',
 model='two-site four fermion modes; interacting 16 Fock states; Ohmic harmonic boson thermal bath; Born-Markov-secular toy; no Lamb shift',
 beta=beta,mu=mu,ohmic_g2=g2,ohmic_cutoff=cutoff,
 oscillator_energy=1.25,oscillator_levels=48,
 number_sector_sizes=[int(sum(ns==i)) for i in range(5)],
 bohr_transition_jump_count=spec['transition_count'],
 collective_zero_frequency_jump_active=bool(np.linalg.norm(spec['zero_frequency_jump'])>1e-6),
 oscillator_correlation_trace_vs_spectral_max=float(abs(corr_mat-corr_pair)),
 oscillator_emission_vs_bose_max=float(abs(os['emission']-(os['occupation']+1))),
 oscillator_absorption_vs_bose_max=float(abs(os['absorption']-os['occupation'])),
 oscillator_KMS_ratio_max=float(abs(os['absorption']/os['emission']-np.exp(-beta*1.25))),
 continuum_boson_KMS_ratio_max=float(max(abs(oscillator_fourier_rate(w,beta=beta,g2=g2,cutoff=cutoff)/oscillator_fourier_rate(-w,beta=beta,g2=g2,cutoff=cutoff)-np.exp(-beta*w)) for w in (.04,.27,1.2,4.5))),
 interacting_Fock_rate_KMS_residual_max=float(np.max(abs(W*p[None,:]-W.T*p[:,None]))),
 interacting_Fock_thermal_stationarity_max=float(np.max(abs(dissipator(rho,jumps)))),
 local_bath_charge_source_abs_max=float(np.max(abs(local['bath_source']))),
 site_charge_continuity_with_bath_source_max=float(np.max(abs(local['continuity']))),
 wrong_omitted_bath_source_local_ward_max=float(np.max(abs(local['density_rate']+np.array([local['current'],-local['current']])))),
 joint_manybody_jump_gauge_covariance_max=float(np.max(abs(dissipator(R(X),Rjs)-R(dissipator(X,jumps))))),
 wrong_stale_jump_gauge_error_max=float(np.max(abs(dissipator(R(X),jumps)-R(dissipator(X,jumps))))),
 bath_contact_finite_difference_max=float(np.max(abs(contact_fd-bath_contact(rho,jumps,G)))),
 magnetic_thermal_stationarity_contact_residual_max=float(np.max(abs(mag['residual']))),
 wrong_omitted_magnetic_contact_static_max=float(np.max(abs(mag['equilibrium_collision_state']))),
 gauge_finite_z_liouville_vs_exact_rotated_density_max=float(np.max(abs(g['dr']-exact_g))),
 gauge_finite_z_charge_density_max=float(np.max(abs(g['density']))),
 gauge_finite_z_complete_physical_current_abs=float(abs(g['current'])),
 gauge_finite_z_charge_Ward_max=float(np.max(abs(g['ward']))),
 wrong_missing_bath_vertex_dynamic_density_matrix_max=float(np.max(abs(gwrong['dr']-exact_g))),
 magnetic_finite_z_density_matrix_abs_max=float(np.max(abs(b['dr']))),
 magnetic_finite_z_charge_Ward_max=float(np.max(abs(b['ward']))),
 wrong_missing_magnetic_bath_vertex_finite_z_response_max=float(np.max(abs(b['dr']-bwrong['dr']))),
 unique_chemical_potential_proved=False,
 ab_initio_phonon_spectrum=False,
 microscopic_electron_self_energy_and_vertex=False,
 source_paper_KME_and_SI_verified=False
)
for key in ('oscillator_correlation_trace_vs_spectral_max',
 'oscillator_emission_vs_bose_max','oscillator_absorption_vs_bose_max',
 'oscillator_KMS_ratio_max','continuum_boson_KMS_ratio_max',
 'interacting_Fock_rate_KMS_residual_max','interacting_Fock_thermal_stationarity_max',
 'site_charge_continuity_with_bath_source_max',
 'joint_manybody_jump_gauge_covariance_max','bath_contact_finite_difference_max',
 'magnetic_thermal_stationarity_contact_residual_max',
 'gauge_finite_z_liouville_vs_exact_rotated_density_max',
 'gauge_finite_z_charge_density_max','gauge_finite_z_complete_physical_current_abs',
 'gauge_finite_z_charge_Ward_max','magnetic_finite_z_charge_Ward_max'):
 assert report[key]<5e-8,(key,report[key])
for key in ('local_bath_charge_source_abs_max','wrong_omitted_bath_source_local_ward_max',
 'wrong_stale_jump_gauge_error_max','wrong_omitted_magnetic_contact_static_max',
 'wrong_missing_bath_vertex_dynamic_density_matrix_max',
 'magnetic_finite_z_density_matrix_abs_max',
 'wrong_missing_magnetic_bath_vertex_finite_z_response_max'):
 assert report[key]>1e-5,(key,report[key])
assert report['collective_zero_frequency_jump_active']
assert report['number_sector_sizes']==[1,4,6,4,1]
reference=ROOT/'results/frozen/round31_bosonic_davies_dynamic_ward.json'
if reference.exists():
    expected=json.loads(reference.read_text())
    assert report['schema']==expected['schema'] and expected['paper_status']=='PARTIAL'
    for key,value in report.items():
        if isinstance(value,float):
            assert abs(value-expected[key])<2e-8,(key,value,expected[key])
print(json.dumps(report,indent=2))
