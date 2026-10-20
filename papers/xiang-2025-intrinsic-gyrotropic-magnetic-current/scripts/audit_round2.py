#!/usr/bin/env python3
"""Second-round audit: TR-related Dirac Hamiltonians, frequency and SI units.

Does NOT fit unknown factors to figure curves or claim complete figure reproduction.
"""
from pathlib import Path
import sys,json,csv
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from igmc_lab.core import (dirac_response,energy_from_frequency,
                            physical_response_si,rashba_analytic)

f=1.0e12
hbar_omega=energy_from_frequency(f,protocol='cyclic_hz')
omega_wrong=energy_from_frequency(f,protocol='angular_rad_s')
g=2.0
N=201
mu_list=np.array([-1.1,-1.,-.85,-.5,-.25,0.,.25,.5,.85,1.,1.1])
data=[]
for mu in mu_list:
    c=dirac_response(N,mu=float(mu),omega_eV=hbar_omega,alpha='x',beta='y')
    d=dirac_response(N,mu=float(mu),omega_eV=hbar_omega,alpha='x',beta='x')
    cr=dirac_response(N,mu=float(mu),omega_eV=hbar_omega,alpha='x',beta='y',mass_sign=-1)
    dr=dirac_response(N,mu=float(mu),omega_eV=hbar_omega,alpha='x',beta='x',mass_sign=-1)
    data.append({'mu_eV':float(mu),'C_si_A_per_m_T':physical_response_si(c['sigma_C'],g_factor=g),
        'C_tr_partner_si':physical_response_si(cr['sigma_C'],g_factor=g),
        'D_si_A_per_m_T':physical_response_si(d['sigma_D'],g_factor=g),
        'D_tr_partner_si':physical_response_si(dr['sigma_D'],g_factor=g),
        'D_wrong_angular_si':physical_response_si(d['sigma_D_per_omega']*omega_wrong,g_factor=g),
        'Kubo_D_full_si':physical_response_si(d['Kubo_imag'],g_factor=g)})

raw=ROOT/'results/raw';fig=ROOT/'results/figures';raw.mkdir(parents=True,exist_ok=True);fig.mkdir(parents=True,exist_ok=True)
with (raw/'round2_si_symmetry.csv').open('w',newline='') as file:
    writer=csv.DictWriter(file,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
ins=dirac_response(N,0,omega_eV=hbar_omega,alpha='x',beta='x')
metal=dirac_response(N,.85,omega_eV=hbar_omega,alpha='x',beta='y')
record={'protocol':'SI_g2_no_fit; source_SM_VI_1THz; TR_partner_beta_and_M_sign_flipped',
  'si_prefactor_A_m-1_T-1_per_angstrom-1':physical_response_si(1.,g_factor=g),
  'frequency_cyclic_hz':f,'photon_energy_cyclic_eV':hbar_omega,
  'photon_energy_if_angular_eV':omega_wrong,
  'rashba_plateau_si_A_m-1_T-1':physical_response_si(rashba_analytic(0.),g_factor=g),
  'dirac_in_gap_D_si_A_m-1_T-1':physical_response_si(ins['sigma_D'],g_factor=g),
  'dirac_in_gap_D_fullKubo_si':physical_response_si(ins['Kubo_imag'],g_factor=g),
  'dirac_metal_mu_0.85_C_si':physical_response_si(metal['sigma_C'],g_factor=g),
  'author_fig2d_visual_estimate_in_gap_si':0.04,
  'author_fig2d_reading_uncertainty_si_approx':0.005,
  'fig2d_visual_estimate_over_computed_in_gap':0.04/physical_response_si(ins['sigma_D'],g_factor=g),
  'author_curve_status':'visually estimated approximate graph reading only, NOT original tabulated data; >2x amplitude discrepancy unresolved',
  'note':'No fitted normalization has been introduced; data are code outputs on independent mesh'}
(raw/'round2_si_audit.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')

x=np.array([r['mu_eV'] for r in data])
figg,ax=plt.subplots(2,2,figsize=(10,7.5),layout='constrained')
a=ax[0,0];a.plot(x,[r['C_si_A_per_m_T'] for r in data],'o-',label='H (C)');a.plot(x,[r['C_tr_partner_si'] for r in data],'x--',label='TR H (C)');a.set(ylabel=r'$-g\mu_B\sigma^C$ [A/(m T)]',title='(a) C is TR-even');a.legend()
a=ax[0,1];a.plot(x,[r['D_si_A_per_m_T'] for r in data],'o-',label='H (D)');a.plot(x,[r['D_tr_partner_si'] for r in data],'x--',label='TR H (D)');a.set(ylabel=r'$-g\mu_B\sigma^D$ [A/(m T)]',title='(b) D is TR-odd');a.legend()
a=ax[1,0];a.plot(x,[r['D_si_A_per_m_T'] for r in data],'o-',label=r'$f=1$ THz, $\hbar 2\pi f$');a.plot(x,[r['D_wrong_angular_si'] for r in data],'x--',label=r'misread as $\omega=10^{12}$ rad/s');a.set(xlabel=r'$\mu_F$ [eV]',ylabel='D [A/(m T)]',title='(c) No fitted factors');a.legend(fontsize=8)
a=ax[1,1];a.axis('off');a.text(.02,.95,'Independent absolute-scale checks',va='top',fontsize=13,weight='bold',transform=a.transAxes)
a.text(.02,.78,f'Rashba Eq.(16), g=2: {record["rashba_plateau_si_A_m-1_T-1"]:.4f} A/(m T)\nPublisher Fig1 plateau: about 4 A/(m T)\n\nDirac D, in-gap: {record["dirac_in_gap_D_si_A_m-1_T-1"]:.5f} A/(m T)\nPublisher Fig2(d): about 0.04 A/(m T)\n\nMagnitude disagreement retained (no fit).',va='top',fontsize=10,transform=a.transAxes)
figg.suptitle('IGMC round 2: time reversal and physical units (computed curves)')
figg.savefig(fig/'round2_TR_and_SI.png',dpi=170);plt.close(figg)
print(json.dumps(record,indent=2))
