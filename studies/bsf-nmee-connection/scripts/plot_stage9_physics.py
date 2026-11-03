"""Physical control variables and measured responses; see STAGE9_FIGURE_QA.md."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bsf_nmee.sliding_tight_binding import lattice_geometry


def main():
    data = json.loads((ROOT/'results/bilayer_stage9_path.json').read_text())
    coupling = json.loads((ROOT/'results/bilayer_stage9_coupling.json').read_text())
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': .7, 'legend.frameon': False})
    fig, axes = plt.subplots(1, 3, figsize=(7.086614, 2.677165))  # 180 x 68 mm
    fig.subplots_adjust(left=.075, right=.93, bottom=.24, top=.83, wspace=.7)
    blue, orange = '#245A81', '#B45C2B'
    source = []
    # One declared representative path; second coupling is a saved robustness check.
    rows = [r for r in data['path'] if r['t_eV'] == .005]
    x = np.array([r['lambda'] for r in rows])
    p = np.array([r['dipole_over_e_d_per_cell'] for r in rows])
    c = np.array([r['C'][0][0] for r in rows])
    axes[0].plot(x, p/p[0], 'o-', ms=3, lw=1.1, color=orange, label=r'$P_z/P_z(\mathrm{AB})$')
    axes[0].plot(x, c/c[0], 's--', ms=3, mfc='white', lw=1.1, color=blue, label=r'$C_{xx}/C_{xx}(\mathrm{AB})$')
    axes[0].set(xlabel=r'AB $\rightarrow$ BA path, $\lambda$', ylabel='Normalized response',
                xticks=[0,.5,1], xticklabels=['AB','midpoint','BA'], ylim=(-1.13,1.13))
    axes[0].legend(loc='upper right', fontsize=6.5)
    source.extend([['a', float(xi), 'Pz/Pz_AB', float(pi)] for xi,pi in zip(x,p/p[0])])
    source.extend([['a', float(xi), 'Cxx/Cxx_AB', float(ci)] for xi,ci in zip(x,c/c[0])])
    rows = sorted([r for r in data['polar_local'] if r['direction'] == [.23,.41]], key=lambda r:r['q'])
    baseline = next(r for r in rows if r['q'] == 0.)
    a, _ = lattice_geometry()
    scale_pm = np.linalg.norm(np.array([.23,.41])@a)*1000
    x = np.array([r['q']*scale_pm for r in rows])
    p = np.array([(r['dipole_over_e_d_per_cell']/baseline['dipole_over_e_d_per_cell']-1)*100 for r in rows])
    n = np.array([r['difference']/abs(baseline['trace'])*100 for r in rows])
    axes[1].plot(x,p,'o-',ms=3,lw=1.1,color=orange,label=r'$\Delta P_z/P_z(\mathrm{AB})$')
    axes[1].plot(x,n,'s--',ms=3,mfc='white',lw=1.1,color=blue,label=r'$(C_{xx}-C_{yy})/|\mathrm{Tr}C_{\rm AB}|$')
    axes[1].set(xlabel='Shear displacement from AB (pm)',ylabel='Relative change (%)',ylim=(-.9,.85))
    axes[1].legend(loc='upper left',fontsize=6.1,handlelength=1.3)
    source.extend([['b',float(xi),'dipole_change_percent',float(pi)] for xi,pi in zip(x,p)])
    source.extend([['b',float(xi),'anisotropy_percent',float(ni)] for xi,ni in zip(x,n)])
    rows = coupling['scan']
    x = np.array([r['parameters']['tunneling_eV']*1000 for r in rows])
    # At t=0 the total spin is zero: its direction is undefined, not zero degrees.
    spin_rows = [r for r in rows if r['parameters']['tunneling_eV'] > 0]
    spin_x = np.array([r['parameters']['tunneling_eV']*1000 for r in spin_rows])
    angle = np.degrees(np.arctan([r['C'][1][0]/r['C'][0][0] for r in spin_rows]))
    p = np.array([r['dipole_over_e_d_per_cell'] for r in rows])*1e6
    axes[2].plot(spin_x,angle,'s-',color=blue,ms=3,lw=1.1,label='Spin tilt')
    right = axes[2].twinx()
    right.spines['right'].set_visible(True)
    right.plot(x,p,'o--',color=orange,mfc='white',ms=3,lw=1.1,label=r'$D_z$')
    axes[2].set(xlabel=r'Interlayer hopping $t_\perp$ (meV)',
                ylabel=r'Spin tilt from $-x$ (deg)',xticks=[0,10,20,30],ylim=(-14,3))
    right.set(ylabel=r'$D_z/(e d_{\rm layer})$ ($10^{-6}$/cell)',ylim=(-.3,2.88))
    axes[2].yaxis.label.set_color(blue);right.yaxis.label.set_color(orange)
    axes[2].tick_params(axis='y',labelcolor=blue);right.tick_params(axis='y',labelcolor=orange)
    axes[2].text(.05,.95,r'Fixed stacking: $\lambda=1/4$',transform=axes[2].transAxes,fontsize=6.5,va='top')
    source.extend([['c',float(xi),'spin_tilt_degrees',float(ci)] for xi,ci in zip(spin_x,angle)])
    source.extend([['c',float(xi),'dipole_in_1e-6',float(pi)] for xi,pi in zip(x,p)])
    for index,(ax,title) in enumerate(zip(axes,['Polar switching','Local shear','Hybridization control'])):
        ax.axhline(0,color='.7',lw=.5,zorder=0)
        ax.text(-.16,1.13,chr(97+index),transform=ax.transAxes,fontweight='bold',fontsize=9)
        ax.set_title(title,fontsize=8,pad=14)
    output = ROOT/'figures/bilayer_stage9_physics'
    fig.savefig(output.with_suffix('.svg'))
    fig.savefig(output.with_suffix('.pdf'))
    fig.savefig(output.with_suffix('.png'), dpi=600)
    plt.close(fig)
    with (ROOT/'results/bilayer_stage9_figure_source.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['panel','x','quantity','value']);writer.writerows(source)


if __name__ == '__main__':
    main()
