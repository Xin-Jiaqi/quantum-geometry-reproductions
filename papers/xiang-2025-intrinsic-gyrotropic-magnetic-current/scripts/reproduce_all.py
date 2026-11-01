#!/usr/bin/env python3
"""Recreate independent model diagnostics for main Figs. 1/2, with frozen parameters.

These are model-level figures in natural units, *not* pixel-perfect or calibrated
published SI amplitudes. Code records numerical provenance and protocol labels.
"""
from __future__ import annotations
import json, csv, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from igmc_lab.core import (dirac_grid,dirac_kernel,dirac_response,rashba_analytic,
                           rashba_integrate,rashba_kernel,kernel_response)
P=json.loads((ROOT/'parameters.json').read_text())
RAW=ROOT/'results/raw';FIG=ROOT/'results/figures'
RAW.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)

def writejson(fname,obj):
    (RAW/fname).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

# Fig.1 2D continuum; e=ℏ=1 kernel, not SI current or arbitrary g factor.
alpha=3.80998212/1.1
x=np.linspace(-.27,.27,601)
y0=alpha*x*x
energy_min=-.6**2/(4*alpha)
mu_r=np.linspace(energy_min-.012,.015,160)
analytic=np.array([rashba_analytic(z) for z in mu_r])
mu_check=[-.015,0.,.009]
numcheck=[rashba_integrate(z,nrad=1200,nphi=96,kmax=.5) for z in mu_check]
X,Y=np.meshgrid(np.linspace(-.23,.23,200),np.linspace(-.23,.23,200),indexing='ij')
kp=rashba_kernel(X,Y,'x','y');kxx=rashba_kernel(X,Y,'x','x')
mask=X*X+Y*Y<.023**2
zxy=np.ma.array(kp['Z01'],mask=mask);zxx=np.ma.array(kxx['Z01'],mask=mask)
fig,axs=plt.subplots(2,2,figsize=(9.4,7.7),layout='constrained')
ax=axs[0,0];ax.plot(x,y0+.6*np.abs(x),label=r'$\epsilon_+$');ax.plot(x,y0-.6*np.abs(x),label=r'$\epsilon_-$');ax.axhline(-.01,c='gray',ls='--',lw=1);ax.set(xlabel=r'$k_x\,({\rm \AA}^{-1})$',ylabel='Energy (eV)',title='(a) Rashba bands');ax.legend()
ax=axs[0,1];ax.plot(mu_r,analytic,label='Eq. (16), e/ℏ removed');ax.plot(mu_check,numcheck,'o',ms=4,label='direct $z=r\\sigma$ quadrature');ax.axvline(0,c='gray',ls=':',lw=.8);ax.set(xlabel=r'$\mu_F$ (eV)',ylabel=r'$\sigma^C_{xy}$ (model $\AA^{-1}$)',title='(b) Intrinsic conduction');ax.legend(fontsize=8)
scale=30
for ax,Z,name in [(axs[1,0],zxy,r'(c) $Z^{-+}_{xy}$'),(axs[1,1],zxx,r'(d) $Z^{-+}_{xx}$')]:
    mesh=ax.pcolormesh(X,Y,Z,cmap='RdBu_r',vmin=-scale,vmax=scale,shading='auto',rasterized=True)
    ax.set_aspect('equal');ax.set(xlabel=r'$k_x\,({\rm \AA}^{-1})$',ylabel=r'$k_y\,({\rm \AA}^{-1})$',title=name)
    fig.colorbar(mesh,ax=ax,shrink=.75)
fig.suptitle('Xiang et al. 2025, Fig. 1: independent model reproduction (not source artwork)')
fig.savefig(FIG/'fig1_rashba_model.png',dpi=175);plt.close(fig)
writejson('fig1_rashba.json',{'protocol':'PRL_Eq15_16_Rashba_analytic_and_interband_polar_quad','units':'e=hbar=1; sigmaC reciprocal angstrom; plotted physical prefactor -gmuB not applied','parameters':P['rashba'],'mu_test_eV':mu_check,'numeric_sigmaC':numcheck,'analytic_sigmaC':[rashba_analytic(z) for z in mu_check],'band_minimum_eV':energy_min,'panels':['a','b','c','d']})

# Fig.2 periodic 2D lattice Dirac; one eigendecomposition for entire mu scan.
N=181
kx,ky=dirac_grid(N)
pairC=dirac_kernel(kx,ky,'x','y')
pairD=dirac_kernel(kx,ky,'x','x')
mu_arr=np.linspace(-1.25,1.25,201)
c_vals=[];d_slope=[]
for mu in mu_arr:
    c_vals.append(kernel_response(pairC,float(mu),omega_eV=.001)['sigma_C'])
    d_slope.append(kernel_response(pairD,float(mu),omega_eV=.001)['sigma_D_per_omega'])
# Band cut at ky=0, plus hotspot and bare metric.
kline=np.linspace(-np.pi,np.pi,601)
dz=.25+2*.1*(2-np.cos(kline)-1)
mag=np.sqrt((.825*np.sin(kline))**2+dz**2)
fig,axs=plt.subplots(2,2,figsize=(9.4,8.0),layout='constrained')
ax=axs[0,0];ax.plot(kline,mag,c='black');ax.plot(kline,-mag,c='black');ax.axhline(0,c='gray',ls='--',lw=.8);ax.set(xlabel=r'$k_x\,({\rm \AA}^{-1})$',ylabel='Energy (eV)',title='(a) Lattice Dirac bands')
ax=axs[0,1];ax.plot(mu_arr,c_vals);ax.axhline(0,c='gray',ls='--',lw=.8);ax.set(xlabel=r'$\mu_F$ (eV)',ylabel=r'$\sigma^C_{xy}$ (natural 2D units)',title='(b) Conduction Hall IGMC')
# K-resolved Q and weighted Q: paper's Fig2(c), no m2-n2 self term.
ax=axs[1,0]; weighted=2*pairD['Q01']/pairD['gap'];lim=np.quantile(np.abs(weighted),.995)
im=ax.pcolormesh(kx,ky,weighted,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto',rasterized=True)
fig.colorbar(im,ax=ax,shrink=.75);ax.set_aspect('equal');ax.set(xlabel=r'$k_x\,(\AA^{-1})$',ylabel=r'$k_y\,(\AA^{-1})$',title=r'(c) $2Q^{-+}_{xx}/\epsilon_{+-}$')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
axin=inset_axes(ax,width='30%',height='30%',loc='upper right')
qbare=pairD['Q01'];sc=np.quantile(np.abs(qbare),.995)
axin.pcolormesh(kx,ky,qbare,cmap='RdBu_r',vmin=-sc,vmax=sc,shading='auto',rasterized=True)
axin.set_xticks([]);axin.set_yticks([]);axin.set_title('bare Q',fontsize=7)
ax=axs[1,1];ax.plot(mu_arr,d_slope);ax.axhline(0,c='gray',ls='--',lw=.8);ax.set(xlabel=r'$\mu_F$ (eV)',ylabel=r'$\sigma^D_{xx}/\Omega$ (natural units)',title=r'(d) Displacement IGMC, $\Omega\rightarrow0$')
fig.suptitle('Xiang et al. 2025, Fig. 2: independent model reproduction (not calibrated SI)')
fig.savefig(FIG/'fig2_dirac_model.png',dpi=175);plt.close(fig)
writejson('fig2_dirac.json',{'protocol':'PRL_Eq17_dirac_low_frequency_geometric_C_D','mesh':[N,N],'parameters':P['dirac'],'temperature_K':0,'mu_eV':mu_arr.tolist(),'sigmaC_xy':c_vals,'sigmaD_xx_per_omega_eV':d_slope,'units':'e=hbar=1; d2k/(2pi)^2; sigmaD computed per photon-energy Omega; paper SI -gmuB conversion not applied','panels':['a','b','c_weighted_plus_bare_inset','d_normalized_by_Omega']})

rows=[]
for name,mu,alpha,beta in [('Dirac_C_insulator',0,'x','y'),('Dirac_D_insulator',0,'x','x'),('Dirac_C_metal',.85,'x','y'),('Dirac_D_metal',.85,'x','x')]:
    for n in [51,101,151,201,301]:
        data=dirac_response(n,mu=mu,omega_eV=.001,alpha=alpha,beta=beta)
        rows.append({'case':name,'mu_eV':mu,'N':n,'sigmaC':data['sigma_C'],'sigmaD_over_omega':data['sigma_D_per_omega'],'Re_Kubo':data['Kubo_real'],'Im_Kubo':data['Kubo_imag']})
with (RAW/'dirac_convergence.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('Wrote Fig1/Fig2 model diagnostics and raw convergence. This is model-only, not absolute SI/pixel acceptance.')
