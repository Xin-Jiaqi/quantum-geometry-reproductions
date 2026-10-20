#!/usr/bin/env python3
"""Round-3 independent response audit; never edits earlier frozen baselines."""
from __future__ import annotations
import csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from igmc_lab.core import (dirac_response,energy_from_frequency,physical_response_si)
from igmc_lab.independent_checks import finite_field_wilson_derivative,orbital_liouville_kernel

params=json.loads((ROOT/'parameters_round3.json').read_text())
w=energy_from_frequency(params['frequency_hz'],protocol=params['frequency_protocol'])
si=physical_response_si(1.,g_factor=params['g_factor'])
conv=[]
for N in params['wilson_meshes']:
    g=dirac_response(N,0.,alpha='x',beta='x',omega_eV=w)
    entry={'N':N,'Q_kernel_dperomega_Ainv_per_eV':g['sigma_D_per_omega'],
           'Kubo_imag_full_Ainv':g['Kubo_imag']}
    for b in params['field_amplitudes_eV']:
        key=f'wilson_dperomega_Ainv_per_eV_b{b:g}'
        entry[key]=finite_field_wilson_derivative(N,b_eV=b)
        entry[f'wilson_relative_error_b{b:g}']=entry[key]/g['sigma_D_per_omega']-1
    conv.append(entry)

liou=[]
for N in params['liouville_meshes']:
    analytic=dirac_response(N,0.,omega_eV=w,alpha='x',beta='x')
    val=orbital_liouville_kernel(N,omega_eV=w)
    liou.append({'N':N,'orbital_liouville_real':float(val.real),
                 'orbital_liouville_imag':float(val.imag),
                 'Kubo_geometric_imag':analytic['Kubo_imag'],
                 'absolute_imag_error':abs(val.imag-analytic['Kubo_imag'])})

source_vis=params['approximate_source_visual_estimates_A_per_m_T']
rows=[]
for mu in np.linspace(-1.22,1.22,81):
    c=dirac_response(params['mu_mesh_N'],mu=float(mu),omega_eV=w,alpha='x',beta='x')
    main=physical_response_si(c['sigma_D'],g_factor=params['g_factor'])
    full=physical_response_si(c['Kubo_imag'],g_factor=params['g_factor'])
    rows.append({'mu_eV':float(mu),'main_eq12_si_A_per_m_T':main,
                 'full_finite_frequency_Kubo_si':full,
                 'SM_eq34_literal_cos_factor2_diagnostic_ONLY':2*main})
raw=ROOT/'results/raw';raw.mkdir(parents=True,exist_ok=True)
with (raw/'round3_displacement_curve.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
with (raw/'round3_wilson_convergence.csv').open('w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(conv[0]));wr.writeheader();wr.writerows(conv)
record={'protocol':'Paper Eq(12) + finite field Wilson loop + orbital Liouville; source SM p6 literal intermediate is DIAGNOSTIC ONLY',
        'frequency_hbaromega_eV':w,'SI_factor_A_per_m_T_per_Ainv':si,
        'parameters':params,'wilson':conv,'liouville':liou,
        'fig2_d_source_visual_estimate':source_vis,
        'main_fig2_d_gap_si':physical_response_si(dirac_response(201,0.,omega_eV=w,alpha='x',beta='x')['sigma_D']),
        'SM_literal_double_gap_si':2*physical_response_si(dirac_response(201,0.,omega_eV=w,alpha='x',beta='x')['sigma_D']),
        'warnings':['Figure point readings are hand-estimated from publisher figure, not original numerical data.',
                    'The SM Eq(34) intermediate normalization is not independently verified and must not be called a physical alternative.',
                    'No fitted normalization or g factor used.']}
(raw/'round3_independent_audit.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
figdir=ROOT/'results/figures';figdir.mkdir(parents=True,exist_ok=True)
fig,(ax,ax2)=plt.subplots(1,2,figsize=(12.8,4.7),layout='constrained')
ax.plot([r['mu_eV'] for r in rows],np.array([r['main_eq12_si_A_per_m_T'] for r in rows])*100,label='Eq.(12) + independent validation')
ax.plot([r['mu_eV'] for r in rows],np.array([r['SM_eq34_literal_cos_factor2_diagnostic_ONLY'] for r in rows])*100,'--',label='SM Eq.(34) intermediate as printed (diagnostic)')
ax.errorbar(source_vis['mu_eV'],np.array(source_vis['displacement_A_per_m_T'])*100,
            yerr=np.array(source_vis['visual_reading_uncertainty_A_per_m_T'])*100,fmt='o',capsize=3,label='Approx. visual paper Fig. 2(d)')
ax.axhline(0,linewidth=.7,color='gray')
ax.set(xlabel=r'Chemical potential $\mu$ (eV)',ylabel=r'$-g\mu_B\sigma^D_{xx}$ ($10^{-2}$ A / m T)',title='Absolute displacement amplitude / source comparison')
ax.legend(fontsize=8)
ns=np.array([x['N'] for x in conv]); errs=np.array([abs(x['wilson_relative_error_b0.0001']) for x in conv]);li=liou[-1]
ax2.semilogy(ns,errs,'o-',label='Wilson-loop finite-field relative error')
ax2.set(xlabel='Uniform k-mesh N (N x N)',ylabel='Relative error |Wilson / Eq.(12) - 1|',title='Independent finite-field convergence')
ax2.text(.05,.12,f"Liouville vs Kubo (N={li['N']}):\nabsolute error {li['absolute_imag_error']:.2e} Å⁻¹",
         transform=ax2.transAxes,fontsize=10,bbox={'facecolor':'white','alpha':.88,'edgecolor':'.8'})
ax2.grid(alpha=.25)
ax2.legend(fontsize=9)
fig.savefig(figdir/'round3_independent_factor_audit.png',dpi=180)
plt.close(fig)
print(json.dumps({'wilson':conv,'liouville':liou,'fig2_d_gap':record['main_fig2_d_gap_si'],'SM_literal_diagnostic_gap':record['SM_literal_double_gap_si']},indent=2))
