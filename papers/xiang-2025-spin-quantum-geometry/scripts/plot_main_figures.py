"""Main Fig1 from the existing paper-model formulas, with source branches visible."""
from pathlib import Path
import sys
import numpy as np
PAPER = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PAPER), str(PAPER.parents[1] / 'scripts')]
from figure_support import setup, export, metadata, plt
from paper_models import fig1a_ndsc, main_fig1b
from paper_eq13 import intrinsic_reduced_eq13


def calculate(nodes=192, radial=72):
    velocity = .6582119569*100  # hbar*v, meV nm, v=10^5 m/s
    tau = .5/.6582119569  # 0.5 ps / (hbar in meV ps)
    mu_a = np.linspace(-80, 80, 401)
    drude_a = np.array([fig1a_ndsc(m, velocity, velocity, .1*velocity, 10., tau, nodes) for m in mu_a])
    intrinsic_a = np.array([intrinsic_reduced_eq13(m, velocity, velocity, .1*velocity, 10., nr=radial, ntheta=nodes)['total'] for m in mu_a])
    mu_b = np.linspace(-400, 400, 801)
    rows = [main_fig1b(m) for m in mu_b]
    keys = ('drude_SM67_natural', 'intrinsic_maintext_natural', 'intrinsic_SM70_natural',
            'intrinsic_literal_SM_integrand_natural', 'independent_density_matrix_static_natural')
    b = {k: np.array([r[k] for r in rows]) for k in keys}
    return dict(mu_a_meV=mu_a, drude_a=drude_a, intrinsic_a=intrinsic_a, mu_b_meV=mu_b, **b)


if __name__ == '__main__':
    setup()
    out = PAPER / 'figures'
    data = calculate()
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.3), layout='constrained')
    ax[0].plot(data['mu_a_meV'], data['drude_a'], label='Drude Eq. (9)')
    ax[0].plot(data['mu_a_meV'], data['intrinsic_a']*100, label='Reduced intrinsic Eq. (13) ×100')
    ax[0].set(xlabel=r'$\mu$ (meV)', ylabel='Natural coefficient (nm/meV)', title='(a) Nonlinear spin current')
    ax[1].plot(data['mu_b_meV']/1000, data['drude_SM67_natural'], label='Drude SM (67)')
    ax[1].plot(data['mu_b_meV']/1000, data['intrinsic_maintext_natural']*1000, label='Printed main intrinsic ×1000')
    ax[1].set(xlabel=r'$\mu$ (eV)', ylabel=r'Coefficient / $(\mu_B e^2)$ (meV$^{-2}$)', title='(b) Nonlinear spin magnetization')
    for a in ax:
        a.axhline(0, color='.6', lw=.7)
        a.legend(fontsize=8)
        a.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    fig.suptitle('Spin Transport — main Fig. 1 formulas\nNatural units; source absolute-axis discrepancy not fitted away', fontsize=12)
    export(fig, out, 'fig1_spin_transport', data)
    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    for key, label in [('intrinsic_maintext_natural','Printed main text'), ('intrinsic_SM70_natural','Printed SM (70)'),
                       ('intrinsic_literal_SM_integrand_natural','Literal SM integral'),
                       ('independent_density_matrix_static_natural','Static density matrix (different protocol)')]:
        ax.plot(data['mu_b_meV']/1000, data[key], label=label)
    ax.set(xlabel=r'$\mu$ (eV)', ylabel=r'Coefficient / $(\mu_B e^2)$ (meV$^{-2}$)', title='Fig. 1(b) intrinsic: source/protocol comparison')
    ax.legend(fontsize=8)
    export(fig, out, 'fig1b_source_comparison', {k:v for k,v in data.items() if k.endswith('_natural') or k=='mu_b_meV'})
    fine = calculate(384, 120)
    metadata(out, {'citation':'10.1103/14mp-263q', 'tau_ps':.5, 'beta_meV':10., 'delta_meV':40.,
             'fig1a_Drude_max_mesh_difference':float(np.max(abs(data['drude_a']-fine['drude_a']))),
             'fig1a_intrinsic_max_mesh_difference':float(np.max(abs(data['intrinsic_a']-fine['intrinsic_a']))),
             'status':'printed model curves recomputed; absolute SI axis not accepted; intrinsic reduced Eq13 only',
             'display_multipliers':{'fig1a_intrinsic':100,'fig1b_intrinsic':1000},
             'multipliers_reason':'same visibility multipliers explicitly used in the paper; not amplitude fits'})
