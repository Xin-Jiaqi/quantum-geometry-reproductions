"""Recompute quantitative main-text panels; no publisher data or fitted factors.

Fig1c is explicitly a frozen-occupation static diagnostic, not a validated
reproduction of the paper's PT-degenerate intrinsic driven-DC coefficient.
"""
from pathlib import Path
import argparse
import sys
import numpy as np

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path[:0] = [str(PAPER), str(ROOT / 'scripts')]
from figure_support import setup, export, metadata, plt
from nmee_lab.clean_core.models import kane_mele
from nmee_lab.clean_core.core import mesh, intrinsic_kernel, fermi, MU_B
from nmee_lab.ptbackend import cumnas as cu


def km_arrays(model, n, mus, components=('xz', 'zx')):
    """Reuse each k kernel across chemical potentials; retain printed/static branches."""
    k = mesh(model.reciprocal, n)
    pref = model.g_factor * MU_B * abs(np.linalg.det(model.reciprocal)) / (2*np.pi)**2 / len(k)
    values = np.zeros((len(components), 2, len(mus)))
    vtop, cbottom = -np.inf, np.inf
    for start in range(0, len(k), 8192):
        sub = k[start:start+8192]
        for j, ag in enumerate(components):
            e, terms, _ = intrinsic_kernel(model, sub, alpha=ag[0], gamma=ag[1], beta=0)
            vtop = max(vtop, float(e[:, 1].max()))
            cbottom = min(cbottom, float(e[:, 2].min()))
            for b, name in enumerate(('printed', 'full_static')):
                for m, mu in enumerate(mus):
                    values[j, b, m] += pref * np.sum(fermi(e, mu, 10.) * terms[name])
    return values, (vtop, cbottom)


def plot_cumnas(n, out):
    p = cu.Model(temperature_K=100., a_nm=1.)
    q = cu.points(n, offset=(.5, .5))
    # With +lambda*sin(kx) and n_y>0 the nodes are on the negative-kx edge.
    vertices = np.array([[-np.pi, np.pi], [0, 0], [0, np.pi], [-np.pi, np.pi]])
    path = np.concatenate([a + (b-a)*np.linspace(0, 1, 181, endpoint=False)[:, None]
                           for a, b in zip(vertices[:-1], vertices[1:])])
    em, ep, _ = cu.energies(path, np.pi/2, p)
    mus = np.linspace(-.2, .2, 161)
    intrinsic = cu.intrinsic_mu_scan(q, 0., mus, p)
    # Actual 2pi sweep: no symmetry-copy of a precomputed quadrant.
    theta = np.linspace(0, 2*np.pi, 73)
    radius = np.linspace(0, .5, 91)
    polar = np.array([cu.extrinsic_mu_scan(q, a, radius, p) for a in theta])
    fig = plt.figure(figsize=(10, 8), layout='constrained')
    ax = fig.add_subplot(221)
    ax.set_axis_off()
    ax.text(.02, .95, '(a) Response schematic (redrawn)', va='top', weight='bold')
    for x, color, label, sign in [(.25, '#C44E52', 'A', 1), (.75, '#4C72B0', 'B', -1)]:
        ax.scatter([x], [.55], s=800, color=color, alpha=.8)
        ax.annotate('', (x, .55+.22*sign), (x, .55), arrowprops={'arrowstyle': '->', 'lw': 2})
        ax.text(x, .3, label, ha='center', weight='bold')
    ax.text(.5, .1, r'$E_x, B_x\ \longrightarrow\ \delta s_z^A=-\delta s_z^B$', ha='center')
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax = fig.add_subplot(222)
    ax.plot(np.arange(len(path))/181, em, color='#4C72B0')
    ax.plot(np.arange(len(path))/181, ep, color='#4C72B0')
    ax.set(xticks=range(4), xticklabels=['M', r'$\Gamma$', 'X', 'M'], ylabel=r'$E/t$',
           title=r'(b) Eq. (5), $\theta=\pi/2$')
    ax = fig.add_subplot(223)
    for j, label in enumerate(('A', 'B')):
        ax.plot(mus, intrinsic[:, j], label=label)
    ax.axhline(0, color='.6', lw=.6)
    ax.set(xlabel=r'$\mu/t$', ylabel=r'$\Gamma\ [(\mathrm{nm\,V\,T})^{-1}]$',
           title='(c) STATIC DIAGNOSTIC — not paper DC')
    ax.legend()
    ax = fig.add_subplot(224, projection='polar')
    ax.grid(False)
    lim = max(abs(polar.min()), abs(polar.max()))
    m = ax.pcolormesh(theta, radius, polar.T, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto', rasterized=True)
    ax.set_title(r'(d) Printed Eq. (2): $\Gamma_{xy;y}$', pad=18)
    ax.set_rlabel_position(135)
    fig.colorbar(m, ax=ax, shrink=.65, label=r'$\Gamma\ [(\mathrm{nm\,V\,T})^{-1}]$')
    fig.suptitle('NMEE / CuMnAs model — independent recalculation\nExplicit model length a = 1 nm; absolute paper amplitude not calibrated', fontsize=12)
    export(fig, out, 'fig1_cumnas', dict(path_k=path, band_lower=em, band_upper=ep,
           mu_intrinsic=mus, static_intrinsic=intrinsic, theta=theta, mu_polar=radius, extrinsic=polar))
    return {'mesh': n, 'static_A_plus_B_max': float(np.max(abs(intrinsic[:, :2].sum(axis=1)))),
            'extrinsic_reversal_max': float(np.max(abs(polar[:37]+polar[36:]))),
            'a_nm': 1., 'temperature_K': 100., 'tau_fs': 50.,
            'fig1c_status': 'static diagnostic, not PT-degenerate paper Eq1 DC reproduction'}


def plot_km(n, out, soc_points):
    model = kane_mele()
    mus = np.linspace(-.3, .3, 161)
    scan, gap = km_arrays(model, n, mus)
    # K and K' along ky=0, matching the pair of Dirac cones in the source.
    kval = 4*np.pi/(3*np.sqrt(3)*.142)
    x = np.linspace(-1.12*kval, 1.12*kval, 601)
    path = np.column_stack((x, np.zeros_like(x)))
    bands = np.linalg.eigvalsh(model.hamiltonian(path))
    soc = np.linspace(-.03, .03, soc_points)
    soc_scan = np.full((2, 2, len(soc)), np.nan)
    soc_gap = np.full((len(soc), 2), np.nan)
    for j, lam in enumerate(soc):
        if abs(lam) < 1.e-12:
            continue  # exact degeneracies: do not use the nondegenerate kernel
        m = kane_mele(intrinsic_soc=float(lam))
        energies = np.linalg.eigvalsh(m.hamiltonian(mesh(m.reciprocal, n)))
        soc_gap[j] = [energies[:, 1].max(), energies[:, 2].min()]
        if soc_gap[j, 1] <= soc_gap[j, 0]:
            raise ValueError('SOC scan has no insulating gap for this mesh/parameter')
        mu = soc_gap[j].mean()
        soc_scan[:, :, j] = km_arrays(m, n, [mu])[0][:, :, 0]
        print(f'SOC {j+1}/{len(soc)}: {lam*1000:.1f} meV', flush=True)
    # Source Fig2(d) is zx;x. Resolve the narrow valley features with two
    # explicitly labelled zooms, rather than a visually empty whole-BZ mesh.
    dq = np.linspace(-.4, .4, 180)  # nm^-1, relative to each valley centre
    xx, yy = np.meshgrid(dq, dq)
    heatmap_k, heatmap = [], []
    for centre in (-kval, kval):
        pts = np.column_stack((xx.ravel()+centre, yy.ravel()))
        e, terms, _ = intrinsic_kernel(model, pts, alpha='z', gamma='x')
        z = model.g_factor*MU_B*np.sum(fermi(e, sum(gap)/2, 10)*terms['printed'], axis=1)
        heatmap_k.append(pts.reshape(*xx.shape, 2))
        heatmap.append(z.reshape(xx.shape))
    heatmap = np.asarray(heatmap)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')
    ax = axes[0, 0]
    ax.plot(x/10, bands*1000, color='.2', lw=.8)
    for g in gap: ax.axhline(g*1000, color='.5', ls='--', lw=.7)
    ax.set(ylim=(-300, 300), xlabel=r'$k_x\ [\mathrm{\AA}^{-1}]$', ylabel='Energy (meV)', title='(a) Eq. (6) bands')
    ax = axes[0, 1]
    for j, label in enumerate(('xz;x', 'zx;x')):
        ax.plot(mus*1000, scan[j, 0], label=label)
        ax.plot(mus*1000, scan[j, 1], ls=':', alpha=.55, color=f'C{j}')
    ax.axvspan(gap[0]*1000, gap[1]*1000, color='.9', zorder=-1)
    ax.set(xlabel=r'$\mu$ (meV)', ylabel=r'$\Gamma\ [(\mathrm{nm\,V\,T})^{-1}]$', title='(b) Solid: printed; dotted: full static')
    ax.legend()
    ax = axes[1, 0]
    for j, label in enumerate(('xz;x', 'zx;x')):
        ax.plot(soc*1000, soc_scan[j, 0], '.-', label=label)
    ax.axhline(0, color='.5', lw=.7)
    ax.set(xlabel=r'$\lambda_{so}$ (meV)', ylabel=r'$\Gamma\ [(\mathrm{nm\,V\,T})^{-1}]$', title='(c) Printed Eq. (1); gap-centre chemical potential')
    ax.legend()
    ax = axes[1, 1]
    ax.set_axis_off()
    ax.set_title('(d) Printed zx;x integrand: valley zooms')
    limit = np.max(abs(heatmap))
    for j, label in enumerate((r'$K$', r"$K'$")):
        sub = ax.inset_axes([.06+.49*j, .23, .39, .59])
        im = sub.pcolormesh(dq/10, dq/10, heatmap[j], cmap='RdBu_r',
                           vmin=-limit, vmax=limit, shading='auto', rasterized=True)
        sub.set(aspect='equal', xlabel=r'$\delta k_x\ [\mathrm{\AA}^{-1}]$', title=label)
        sub.tick_params(labelsize=8)
        if j == 0:
            sub.set_ylabel(r'$k_y\ [\mathrm{\AA}^{-1}]$')
        else:
            sub.set_yticklabels([])
    fig.colorbar(im, ax=ax, shrink=.7, label='Unintegrated kernel (model units)')
    fig.suptitle('NMEE / Kane–Mele — recomputed formulas, no fitted amplitude\nPublished and full-static protocols kept separate', fontsize=12)
    export(fig, out, 'fig2_kane_mele', dict(mu=mus, response=scan, lambda_so=soc, soc_response=soc_scan,
           soc_gap=soc_gap, soc_mu=soc_gap.mean(axis=1), path_k=path, bands=bands,
           heatmap_k=np.asarray(heatmap_k), heatmap_dk=dq, heatmap=heatmap))
    return {'mesh': n, 'soc_points': soc_points, 'soc_zero': 'not evaluated: nondegenerate formula invalid',
            'gap_eV': list(gap), 'temperature_K': 10., 'bond_nm': .142,
            'soc_mu_protocol': 'mu follows each SOC-dependent gap centre; not a literal fixed-mu source reconstruction',
            'heatmap': 'zx;x printed kernel; two K/K-prime patches, relative k +/-0.04 angstrom^-1',
            'fig2_status': 'printed Eq1 model recomputation; not digitized-curve validation'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, default=451)
    parser.add_argument('--cumnas-mesh', type=int, default=201)
    parser.add_argument('--soc-points', type=int, default=41)
    args = parser.parse_args()
    setup()
    out = PAPER / 'figures'
    checks = {'citation': '10.1103/r3x2-kb93', 'generated_by': 'scripts/plot_main_figures.py'}
    checks['fig1'] = plot_cumnas(args.cumnas_mesh, out)
    checks['fig2'] = plot_km(args.mesh, out, args.soc_points)
    checks['convergence'] = {}
    for n in (args.mesh, args.mesh+150):
        checks['convergence'][str(n)] = km_arrays(kane_mele(), n, [.025])[0][:, 0, 0].tolist()
    metadata(out, checks)
