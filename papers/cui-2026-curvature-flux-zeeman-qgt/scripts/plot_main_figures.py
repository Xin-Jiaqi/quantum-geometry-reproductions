"""Redraw the two conceptual main figures and add a computed flux check.

The paper's Figs1/2 are schematics, not numerical transport curves.
The inward OmegaA arrows follow its signed Eq8; -2 OmegaA would point outward.
"""
from pathlib import Path
import sys
import numpy as np
PAPER=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(PAPER/'src'),str(PAPER.parents[1]/'scripts')]
from figure_support import setup,export,metadata,plt
from flux_lab.core import model,qgt_eigen,qgt_projector,winding_field,flux_circle_from_qgt,berry_phase_wilson_circle


if __name__=='__main__':
    setup()
    out=PAPER/'figures'
    fig,ax=plt.subplots(figsize=(10,4),layout='constrained')
    ax.set_axis_off()
    ax.text(.5,.93,'Zeeman QGT: $T^{ab}=r^a\\sigma^b$',ha='center',fontsize=15)
    labels=[('Normal metric',r'$g^N$: real symmetric'),('Normal curvature',r'$\Omega^N$: imaginary antisymmetric'),
            ('Anomalous metric',r'$g^A$: imaginary symmetric'),('Anomalous curvature',r'$\Omega^A$: real antisymmetric')]
    for x,(name,label) in zip(np.linspace(.12,.88,4),labels):
        ax.annotate('',(x,.63),(.5,.84),arrowprops={'arrowstyle':'->','color':'.5'})
        ax.text(x,.50,name,ha='center',fontsize=10,weight='bold')
        ax.text(x,.35,label.replace(': ',':\n'),ha='center',fontsize=9)
    ax.text(.5,.08,'Fig. 1 redrawn: symmetry under index exchange defines metric vs curvature',ha='center')
    export(fig,out,'fig1_geometry_decomposition',{'sector_order':np.array([x[0] for x in labels])})

    theta=np.arange(16)*2*np.pi/16
    xy=np.column_stack((np.cos(theta),np.sin(theta)))
    omega=np.array([qgt_eigen(*model(*k))['omegaA'][:2] for k in xy])
    winding=np.array([winding_field(*k) for k in xy])
    fig,ax=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    ax[0].add_patch(plt.Circle((0,0),1,fill=False,color='.7'))
    ax[0].quiver(xy[:,0],xy[:,1],winding[:,0],winding[:,1],color='#4C72B0',angles='xy',scale_units='xy',scale=3,label='w: winding field')
    ax[0].quiver(xy[:,0],xy[:,1],omega[:,0],omega[:,1],color='#DD8452',angles='xy',scale_units='xy',scale=1.5,label=r'$\Omega^A$: flux field')
    ax[0].scatter([0],[0],s=35,color='.2')
    ax[0].set(aspect='equal',xlim=(-1.5,1.5),ylim=(-1.5,1.5),xlabel=r'$k_x$',ylabel=r'$k_y$',title='(a) QGT-computed vectors; signed Eq. (8)')
    ax[0].legend(loc='lower center',fontsize=8)
    ax[1].set_axis_off()
    table=ax[1].table(cellText=[[r'Im $\sigma^-$',r'$\omega$',r'$\langle\Omega^A\rangle$'],
                               [r'Im $\sigma^+$',r'$\omega$',r'$\langle g^N\rangle$'],
                               [r'Re $\sigma^-$',r'$\omega^2$',r'$\langle\Omega^N\rangle$'],
                               [r'Re $\sigma^+$',r'$\omega^2$',r'$\langle g^A\rangle$']],
                        colLabels=['Channel','Scaling','Geometry'],cellLoc='center',loc='center')
    table.scale(1,2.5)
    ax[1].set_title('(b) Source Eq. (18): dynamic, static-subtracted response')
    ax[1].text(.5,.08,'A theoretical channel map; not simulated measurements',ha='center',fontsize=9)
    fig.suptitle('Curvature–Flux — main Fig. 2 redrawn and vector field recomputed',fontsize=12)
    export(fig,out,'fig2_flux_and_channels',dict(k=xy,winding=winding,omegaA=omega))

    radii=np.geomspace(.05,3,25)
    values=np.array([flux_circle_from_qgt(r,points=512) for r in radii])
    fine=np.array([flux_circle_from_qgt(r,points=1024) for r in radii])
    massive=np.array([flux_circle_from_qgt(r,m=.25,points=512) for r in radii])
    analytic=-np.pi*radii**2/(radii**2+.25**2)
    wilson=np.array([abs(berry_phase_wilson_circle(r,points=512))/np.pi for r in radii])
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    ax.semilogx(radii,-values/np.pi,'o-',label=r'Gapless: $-\Phi/\pi$ from QGT')
    ax.semilogx(radii,-massive/np.pi,'s-',label=r'Gapped $m=0.25$: QGT flux')
    ax.semilogx(radii,-analytic/np.pi,':',color='.3',label='Gapped analytic expression')
    ax.semilogx(radii,wilson,'x',label=r'Gapless Wilson loop $|\gamma|/\pi$')
    ax.set(xlabel='Contour radius (model units)',ylabel='Normalized invariant / flux',title='Additional numerical check of the main-text topology claim')
    ax.legend(fontsize=8)
    export(fig,out,'flux_numerical_check',dict(radius=radii,flux=values,flux_refined=fine,massive_flux=massive,massive_analytic=analytic,wilson_abs_over_pi=wilson))
    d,dk=model(.3,.4,.2)
    error=np.max(abs(qgt_eigen(d,dk)['T']-qgt_projector(d,dk)['T']))
    metadata(out,{'citation':'arXiv:2604.09725v3','source_figures':'Figs1/2 are conceptual schematics',
                  'added_check':'eigenvector-QGT contour flux; independent Wilson-loop topology',
                  'matrix_vs_projector_error':float(error),'mesh_refinement_max_error':float(max(abs(values-fine))),
                  'gapless_flux_error_vs_minus_pi':float(max(abs(fine+np.pi))),
                  'massive_flux_error':float(max(abs(massive-analytic))),
                  'sign':'OmegaA inward for d=(kx,ky,0), flux=-pi. Outward vector would be -2 OmegaA.',
                  'scope':'local model topology and source channel table; not a bulk current or experimental curve'})
