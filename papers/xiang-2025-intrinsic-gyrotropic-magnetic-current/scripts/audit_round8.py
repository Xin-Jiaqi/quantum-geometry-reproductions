#!/usr/bin/env python3
"""Round 8: SM semiclassical/Green's function independent closure.

With --no-write prints JSON only; otherwise writes raw results & one figure.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from igmc_lab.core import dirac_grid, dirac_hamiltonian,energy_from_frequency,pair_kernels,physical_response_si
from igmc_lab.green_wavepacket_audit import (green_eq44_resolvent,wavepacket_cos_sin,
                                              analytical_cosine_harmonics)

p=argparse.ArgumentParser();p.add_argument('--no-write',action='store_true');a=p.parse_args()
parameters=json.loads((ROOT/'parameters_round8.json').read_text())
n=parameters['mesh'];kx,ky=dirac_grid(n)
H,dH=dirac_hamiltonian(kx,ky)
w=energy_from_frequency(parameters['frequency_hz'],protocol='cyclic_hz')
results={'response_protocol':'SM_Eq21_24_wavepacket_AND_SM_Eq44_lesserGreen_resolvent',
         'physical_scope':'two-band periodic Dirac, filled lower band, no disorder, nonresonant',
         'mesh':[n,n], 'photon_energy_eV':w,'eta_eV':0.,'chemical_potential_eV':0.,
         'columns':{},'SM_eq46_FT_coefficient_positive_harmonic':analytical_cosine_harmonics()[0],
         'SM_eq46_printed_over_correct':2., 'SM_eq34_metric_intermediate_over_main':2.}
max_res=0.;max_stat=0.; max_q=0.; max_wrong_sine=0.
for axis,spin in [('x','x'),('x','y'),('y','x'),('y','y'),('x','z'),('y','z')]:
    wave=wavepacket_cos_sin(H,dH[axis],spin,w)
    green=green_eq44_resolvent(H,dH[axis],spin,w)
    delta=np.max(abs(green-wave['K_from_wavepacket']))
    max_res=max(max_res,float(delta))
    # Direct model spectral kernel only as an additional control; not used by either new method.
    pair=pair_kernels(H,dH,axis,spin)
    baseline=pair['amplitude01']/(w-pair['gap']) -pair['amplitude01'].conj()/(w+pair['gap'])
    max_stat=max(max_stat,float(np.max(abs(green-baseline))))
    results['columns'][f'{axis}{spin}']={
       'green_mean_K_Ainv':{'real':float(np.mean(green.real)), 'imag':float(np.imag(np.mean(green)))},
       'wave_mean_K_Ainv':{'real':float(np.mean(wave['K_from_wavepacket'].real)), 'imag':float(np.mean(wave['K_from_wavepacket'].imag))},
       'max_abs_green_vs_wave_per_k':float(delta),
       'P_cos_over_b_Ainv':float(np.mean(wave['P_cos_per_b'])),
       'P_sin_over_b_Ainv':float(np.mean(wave['P_sin_per_b'])),
       'P_sin_main_Eq10_over_b_Ainv':float(np.mean(wave['P_sin_main_Eq10_per_b'])),
       'physical_imag_K_SI_A_m_inv_T_inv':physical_response_si(float(np.imag(np.mean(green))))
    }
    max_wrong_sine=max(max_wrong_sine,float(np.max(abs(wave['P_sin_per_b']-wave['P_sin_main_Eq10_per_b']))))
results['max_per_k_GF_vs_wave_abs']=max_res
results['max_per_k_GF_vs_spectral_abs']=max_stat
results['max_per_k_mainEq10_F_sine_difference']=max_wrong_sine
results['fig2d_Dxx_1THz_GF_SI']=results['columns']['xx']['physical_imag_K_SI_A_m_inv_T_inv']
results['source_notes']={
 'SM_eq34':'SM page 6 Eq.(34) writes Q=2Re(r sigma) after starting with 2Re of the wavepacket and B harmonics each carrying 1/2. The correct coefficient for its displayed denominator bracket is Re(r sigma).',
 'SM_eq34_F_vs_main10':'Main Eq.(10) says F=-omega Z/gap^2 and P=-gmuB F sin; direct SM Eq.(21,24) gives P_sin/b=-omega Z/(gap^2-omega^2). Their signs disagree at O(omega); its derivative first affects J at O(omega^2), omitted in Eq(12).',
 'SM_eq46':'SM page 7 Eq.(46) writes the Fourier integral of cos with 2pi delta at each +/- frequency; with declared integral convention coefficient is pi. Eq.(47) actually uses the correct harmonic field amplitude B/2.',
 'paper_fig2d_status':'UNRESOLVED original numeric data not supplied; do not infer author plot from local factor-2 algebra alone',
 'main_eq14_status':'Open response-domain antiunitary convention; round6 already tested K_TR(omega)=K_H(omega)* in the real phase convention.'
}
if not a.no_write:
    folder=ROOT/'results/raw';folder.mkdir(parents=True,exist_ok=True)
    (folder/'round8_GF_wavepacket_audit.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # One generic k-point where all quadratures are nonzero.
    h,dh=dirac_hamiltonian(np.array([.47]),np.array([-.58]))
    freqs=np.linspace(.001,.12,80)
    pcos=[];psin=[];sin_main=[];diff_GF=[]
    for ww in freqs:
        wave=wavepacket_cos_sin(h,dh['x'],'x',ww)
        pcos.append(float(wave['P_cos_per_b'][0]))
        psin.append(float(wave['P_sin_per_b'][0]))
        sin_main.append(float(wave['P_sin_main_Eq10_per_b'][0]))
        gf=green_eq44_resolvent(h,dh['x'],'x',ww)
        diff_GF.append(float(abs(gf[0]-wave['K_from_wavepacket'][0])))
    fig,axs=plt.subplots(2,1,figsize=(8.2,6.6),layout='constrained')
    axs[0].plot(freqs,psin,label='SM Eq.(21,24), exact P_sin/b')
    axs[0].plot(freqs,sin_main,ls='--',label='Main Eq.(10) F sinusoid')
    axs[0].set(xlabel=r'Photon energy $\hbar\omega$ (eV)',ylabel=r'$P_{\sin}/b$ (model units)',title='Wavepacket sine-term sign check at one k point')
    axs[0].legend(fontsize=9)
    axs[1].semilogy(freqs,np.maximum(diff_GF,1e-17),label='|GF Eq.44 - wavepacket Eqs.21,24|')
    axs[1].set(xlabel=r'Photon energy $\hbar\omega$ (eV)',ylabel='Absolute difference',title='Independent finite-frequency response equality')
    axs[1].legend(fontsize=9)
    fig.savefig(ROOT/'results/figures/round8_GF_wavepacket_audit.png',dpi=170);plt.close(fig)
print(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False))
