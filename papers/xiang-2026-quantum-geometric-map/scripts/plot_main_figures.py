"""Recompute main Fig2 with literal Eq11 amplitude and documented finite-T continuation."""
from pathlib import Path
import sys
import numpy as np
PAPER = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PAPER/'src'), str(PAPER.parents[1]/'scripts')]
from figure_support import setup, export, metadata, plt
from qgmap.dirac import Parameters, energies, local_zbc_analytic, local_zbc_eigensolver, continuum_temperature_kernel, zeeman_time_prefactor
from qgmap.fig2_si_v52 import absolute_fig2_audit


if __name__ == '__main__':
    setup()
    out = PAPER/'figures'
    p = Parameters()
    x = np.linspace(-.06, .06, 501)
    em, ep = energies(x, 0., p)
    mu = np.linspace(-.04, .04, 801)
    curves = []
    for temperature in (0, 20, 40):
        y = np.full_like(mu, np.nan)
        valid = mu != 0 if temperature==0 else np.ones(len(mu), bool)
        y[valid] = -1e4*zeeman_time_prefactor(p)*continuum_temperature_kernel(mu[valid],temperature,p.beta)
        curves.append(y)
    grid = np.linspace(-.05,.05,242)  # node excluded, no division by zero
    xx, yy = np.meshgrid(grid,grid)
    zx, zy = local_zbc_analytic(xx, yy)
    fig, ax = plt.subplots(2,2,figsize=(10,8),layout='constrained')
    ax[0,0].plot(x,em*1000,color='.2');ax[0,0].plot(x,ep*1000,color='.2')
    ax[0,0].axhline(10,color='.5',ls='--')
    ax[0,0].set(ylim=(-40,40),xlabel=r'$k_x$ ($\mathrm{\AA}^{-1}$)',ylabel='Energy (meV)',title='(a) Tilted Dirac bands, Eq. (10)')
    for T,y in zip((0,20,40),curves):ax[0,1].plot(mu*1000,y,label=f'{T} K')
    ax[0,1].axhline(1.75,color='.5',ls=':',label='Published plateau ≈1.75')
    ax[0,1].set(xlabel=r'$\mu$ (meV)',ylabel=r'$-10^4\Delta\sigma/(e^2/\hbar)$',title='(b) Literal Eq. (11); finite-T continuation')
    ax[0,1].legend(fontsize=8)
    for a,z,label in zip(ax[1],(zx,zy),(r'(c) $Z_{-+}^{xx}$',r'(d) $Z_{-+}^{yx}$')):
        # Symmetric log keeps every finite point visible without clipping the node hotspot.
        from matplotlib.colors import SymLogNorm
        lim=float(np.max(abs(z)))
        im=a.pcolormesh(grid,grid,z,cmap='RdBu_r',norm=SymLogNorm(linthresh=20,vmin=-lim,vmax=lim),shading='auto',rasterized=True)
        a.set(aspect='equal',xlabel=r'$k_x$ ($\mathrm{\AA}^{-1}$)',ylabel=r'$k_y$ ($\mathrm{\AA}^{-1}$)',title=label)
        fig.colorbar(im,ax=a,shrink=.7,label='Zeeman curvature (Å); symlog color')
    fig.suptitle('Magnetotransport — main Fig. 2 recalculation\nBand/geometry model; response amplitude discrepancy shown explicitly',fontsize=12)
    export(fig,out,'fig2_magnetotransport',dict(k=x,band_lower=em,band_upper=ep,mu=mu,temperatures=np.array([0,20,40]),response=curves,grid=grid,zxx=zx,zyx=zy))
    # Eq11 amplitude checked by an independent SI restoration already in repository.
    audit=absolute_fig2_audit()
    checks=[]
    for kx,ky in ((.01,.02),(-.02,.03),(.04,-.01)):
        checks.append(float(np.max(abs(np.asarray(local_zbc_analytic(kx,ky))-local_zbc_eigensolver(kx,ky)))))
    metadata(out,{'citation':'10.1103/kd12-4ycj','parameters':p.__dict__,
                  'geometry_eigenvector_max_error':max(checks),'SI_audit':audit,
                  'status':'Fig2 numerical panels recomputed; 2pi amplitude mismatch remains explicit',
                  'finite_temperature':'tanh Fermi-population continuation; printed Eq11 is T=0',
                  'fig1':'conceptual classification diagram, not numerical data; summarized in gallery'})
