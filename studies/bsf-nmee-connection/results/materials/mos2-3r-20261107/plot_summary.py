"""Plot all seven planned MoS2 control points, without fitting a switch field."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / 'summary.json').read_text())
rows = {(r['structure'], r['gate_Vnm']): r for r in data['rows']
        if r['mesh'] == 6 and r['probe_VA'] == .004 and r['zeeman_eV'] == .016}
# The smaller-probe check is reported separately, not mixed into the gate curve.
required = {('mono', f) for f in (-.25, 0., .25)} | {('A', f) for f in (-.25, 0., .25)} | {('B', 0.)}
if not required <= rows.keys():
    raise SystemExit(f'Incomplete planned comparison: {required - rows.keys()}')
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .7, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
})
fig, axes = plt.subplots(1, 2, figsize=(7.007874, 3.070866))  # 178 x 78 mm
fig.subplots_adjust(left=.10, right=.98, bottom=.22, top=.82, wspace=.38)
colors = ['#3979A8', '#D88A35', '#252525']
markers = ['o', 's', 'D']
fields = ['C_bottom', 'C_top', 'C_total']
labels = ['Lower layer', 'Upper layer', 'Net']
unit = r'$C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]'
for ax, letter, title in zip(axes, 'ab', ['Zero perpendicular bias', 'Fixed 3R stacking A']):
    ax.axhline(0, color='.7', lw=.7, zorder=0)
    ax.set_ylabel(unit)
    ax.text(-.18, 1.10, letter, fontweight='bold', transform=ax.transAxes)
    ax.set_title(title, fontsize=9, pad=12)
ax = axes[0]
ax.scatter(0, rows['mono', 0.]['C_total'], color=colors[2], marker='D', s=30, zorder=4)
for x, name in [(1, 'A'), (2, 'B')]:
    r = rows[name, 0.]
    for dx, key, color, marker in zip([-.16, 0, .16], fields, colors, markers):
        ax.scatter(x + dx, r[key], color=color, marker=marker, s=30, zorder=4)
ax.set_xticks([0, 1, 2], ['Isolated\nmonolayer', '3R A', 'Mirror 3R B'])
ax.set_xlim(-.45, 2.5)
ax = axes[1]
gates = np.array([-.25, 0., .25])
for key, color, marker, label in zip(fields, colors, markers, labels):
    ax.plot(gates, [rows['A', f][key] for f in gates], color=color,
            marker=marker, ms=4, lw=1.2, label=label)
ax.plot(gates, [2 * rows['mono', f]['C_total'] for f in gates],
        color='#888888', marker='^', ms=4, ls='--', lw=1,
        label='2 × isolated monolayer')
ax.set_xlabel(r'Perpendicular bias $F_z$ (V/nm)')
ax.set_xticks(gates)
ax.legend(loc='best', frameon=False)
fig.savefig(ROOT / 'MoS2_NMEE_comparison.svg')
fig.savefig(ROOT / 'MoS2_NMEE_comparison.pdf')
fig.savefig(ROOT / 'MoS2_NMEE_comparison.png', dpi=600)
plt.close(fig)
