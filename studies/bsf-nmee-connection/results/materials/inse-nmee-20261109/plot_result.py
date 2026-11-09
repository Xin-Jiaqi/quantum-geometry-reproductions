"""Interface electrical gain and its signed layer distribution; deterministic data."""
from pathlib import Path
import json, csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

p = Path(__file__).resolve().parent
s = json.loads((p/'summary.json').read_text()); assert s['checks']['all_seven_complete']
ref = json.loads((p/'reference.json').read_text()); assert ref['complete']
mo = json.loads((p/'mosse_comparison.json').read_text())['slopes']
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 8,
    'axes.titlesize': 8.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': .65, 'lines.linewidth': 1.2,
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'savefig.dpi': 600})
blue, orange, green, gray, ink = '#0072B2', '#D55E00', '#009E73', '#777777', '#222222'
fig, axs = plt.subplots(1, 3, figsize=(7.204724409, 2.637795276))
fig.subplots_adjust(left=.10, right=.985, bottom=.24, top=.84, wspace=.55)
fields = np.array([-.1, 0., .1]); source = []

def point(name, nb, key):
    return next(r[key] for r in s['points'] if r['name'] == name and r['nb'] == nb)

def actual_curve(stem, nb, key='C_net', factor=1.):
    names = ['Am', 'A0', 'Ap'] if stem == 'A' else ['monom', 'mono0', 'monop']
    c = np.array([point(n, nb, key) for n in names])*factor
    return c-c[1]

def reference_curve(nb, end):
    vals = []
    for case in ('Am', 'A0', 'Ap'):
        vals.append(sum(r['C_layer'] for r in ref['points'][case]['rows']
                        if r['nb'] == nb and round(r['taper_end_A'], 6) == end))
    vals = np.array(vals); return vals-vals[1]

def curves(ax, panel, label, values, color, marker, selected=-1, linestyle='-'):
    values = np.array(values)
    ax.fill_between(fields, values.min(axis=0), values.max(axis=0), color=color, alpha=.13, linewidth=0)
    ax.plot(fields, values[selected], marker=marker, markersize=3.8, color=color, linestyle=linestyle, label=label)
    for j, y in enumerate(values):
        source.extend(dict(panel=panel, series=label, setting=j, x=float(x), y=float(v), representative=j == selected % len(values)) for x, v in zip(fields, y))

curves(axs[0], 'a', 'Bilayer', [actual_curve('A', n) for n in (128, 144, 160)], ink, 'o')
curves(axs[0], 'a', 'Isolated pair', [actual_curve('mono', n, factor=2) for n in (64, 72, 80)], gray, 's', linestyle='--')
ends = sorted({round(r['taper_end_A'], 6) for r in ref['points']['A0']['rows']})
curves(axs[0], 'a', 'Electrostatic', [reference_curve(n, e) for n in (64, 72, 80) for e in ends], green, '^', linestyle=':')
axs[0].set_title('Net electrical control', pad=9)
axs[0].set_ylabel(r'$\Delta C_{\rm net}$'+'\n'+r'[$10^{-9}$ (nm V T)$^{-1}$]')
axs[0].legend(frameon=False, loc='upper right', handlelength=1.6, labelspacing=.28)

curves(axs[1], 'b', 'Lower layer', [actual_curve('A', n, 'C_lower') for n in (128, 144, 160)], orange, 'v')
curves(axs[1], 'b', 'Upper layer', [actual_curve('A', n, 'C_upper') for n in (128, 144, 160)], blue, '^')
curves(axs[1], 'b', 'Isolated layer', [actual_curve('mono', n) for n in (64, 72, 80)], gray, 's', linestyle='--')
axs[1].set_title('Changes in each layer', pad=9)
axs[1].set_ylabel(r'$\Delta C_{\ell}$'+'\n'+r'[$10^{-9}$ (nm V T)$^{-1}$]')
axs[1].legend(frameon=False, loc='upper right', handlelength=1.6, labelspacing=.28)
for ax in axs[:2]:
    ax.axhline(0, color='#BDBDBD', linewidth=.55, zorder=0)
    ax.set_xlabel(r'Gate field $F$ (V nm$^{-1}$)')
    ax.set_xticks(fields); ax.set_xlim(-.112, .112)

inse_g = next(r for r in s['slopes'] if r['system'] == 'A' and r['nb'] == 160)
inse_pair = next(r['matched_pair_G'] for r in s['slopes'] if r['system'] == 'mono' and r['nb'] == 80)
mo_pair = mo['mono']['64'][0]+mo['upper']['64'][0]
ratios = np.array([[.5, .5], np.array(mo['A']['112'][1:])/mo_pair,
                   np.array([inse_g['G_lower'], inse_g['G_upper']])/inse_pair])
assert abs(inse_pair) > 1e-4 and abs(mo_pair) > 1e-4
ax = axs[2]
for i, pair in enumerate(ratios):
    pos = neg = 0.
    for name, value, color in zip(('Lower layer', 'Upper layer'), pair, (orange, blue)):
        bottom = pos if value > 0 else neg
        ax.bar(i, value, bottom=bottom, width=.57, color=color, edgecolor='white', linewidth=.5)
        if value > 0: pos += value
        else: neg += value
        source.append(dict(panel='c', series=name, setting='representative', x=i, y=float(value), representative=True))
    ax.plot(i, sum(pair), 'D', color=ink, markeredgecolor='white', markeredgewidth=.6, markersize=4.8, label='Net gain' if i == 0 else None)
    source.append(dict(panel='c', series='Net gain', setting='representative', x=i, y=float(sum(pair)), representative=True))
    ax.annotate(f'{sum(pair):.2f}', (i, sum(pair)), xytext=(8, 0), textcoords='offset points', va='center', fontsize=7)
alternate = {
    1: [v[0]/(mo['mono'][nm][0]+mo['upper'][nm][0]) for v in mo['A'].values() for nm in ('48', '56', '64')],
    2: [a['G_net']/m['matched_pair_G'] for a in s['slopes'] if a['system'] == 'A'
        for m in s['slopes'] if m['system'] == 'mono']}
for i, vals in alternate.items():
    ax.plot(i-.2+np.linspace(-.025, .025, len(vals)), vals, 'o', ms=2.4, mfc='white', mec=gray, mew=.5, linestyle='none')
    source.extend(dict(panel='c', series='Net gain', setting=j, x=i, y=float(v), representative=False) for j, v in enumerate(vals))
ax.axhline(0, color='#BDBDBD', lw=.55)
ax.set_xticks([0, 1, 2], ['Isolated\npair', 'MoSSe', 'InSe'])
ax.set_xlim(-.55, 2.75); ax.set_ylim(-3.25, 4.8)
ax.set_ylabel(r'$G_{\ell}/G_{\rm isolated\ pair}$')
ax.set_title('Different layer contributions', pad=9)
ax.legend(frameon=False, loc='lower left', handlelength=1.0)
for label, ax in zip('abc', axs):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(direction='out', length=2.5, width=.6)
    ax.text(-.22, 1.14, label, transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
fig.savefig(p/'inse_gate_comparison.png', dpi=600, facecolor='white')
fig.savefig(p/'inse_gate_comparison.svg', facecolor='white')
fig.savefig(p/'inse_gate_comparison.pdf', facecolor='white')
svg = p/'inse_gate_comparison.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
with (p/'figure_source.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(source[0]), lineterminator='\n'); w.writeheader(); w.writerows(source)
plt.close(fig)
