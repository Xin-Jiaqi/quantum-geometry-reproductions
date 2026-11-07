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
fig, axes = plt.subplots(1, 3, figsize=(7.007874, 3.070866))  # 178 x 78 mm
fig.subplots_adjust(left=.09, right=.985, bottom=.24, top=.75, wspace=.52)
colors = ['#3979A8', '#D88A35', '#252525']
markers = ['o', 's', 'D']
fields = ['C_bottom', 'C_top', 'C_total']
labels = ['Lower layer', 'Upper layer', 'Net']
unit = r'Response coefficient [$10^{-9}$ (nm V T)$^{-1}$]'
for ax, letter, title in zip(axes, 'abc', ['Interface activation', r'Gate-induced change $\Delta C_{xx}$', 'Net response']):
    ax.axhline(0, color='.7', lw=.7, zorder=0)
    ax.text(-.18, 1.10, letter, fontweight='bold', transform=ax.transAxes)
    ax.set_title(title, fontsize=8, pad=12)
fig.text(.008, .50, unit, va='center', rotation=90, fontsize=8)
ax = axes[0]
ax.scatter(0, rows['mono', 0.]['C_total'], color=colors[2], marker='D', s=30, zorder=4)
for x, name in [(1, 'A'), (2, 'B')]:
    r = rows[name, 0.]
    for dx, key, color, marker in zip([-.16, 0, .16], fields, colors, markers):
        ax.scatter(x + dx, r[key], color=color, marker=marker, s=24, zorder=4,
                   label=labels[fields.index(key)] if x==1 else None)
ax.set_xticks([0, 1, 2], ['Isolated\nmonolayer', '3R A', '3R B'])
ax.set_xlim(-.45, 2.5)
ax = axes[1]
gates = np.array([-.25, 0., .25])
for key, color, marker, label in zip(fields[:2], colors[:2], markers[:2], labels[:2]):
    ax.plot(gates, [rows['A', f][key]-rows['A',0.][key] for f in gates], color=color,
            marker=marker, ms=4, lw=1.2, label=label)
ax.set_xlabel(r'Perpendicular bias $F_z$ (V/nm)')
ax.set_xticks(gates)
ax=axes[2]
ax.plot(gates, [rows['A', f]['C_total'] for f in gates],
        color=colors[2], marker='D', ms=4, lw=1.2, label='Bilayer A')
forecast=json.loads((ROOT/'gate_prediction.json').read_text())
pred={r['gate_physical_Vnm']:r['predicted_C_total'] for r in forecast['predictions']}
pred[0.]=rows['A',0.]['C_total']
ax.plot(gates, [pred[f] for f in gates],color='#89769A',marker='v',ms=4,
        mfc='white',ls='--',lw=1,label='Additive prediction')
ax.plot(gates, [2*rows['mono',f]['C_total'] for f in gates],
        color='#888888', marker='^', ms=4, ls=':', lw=1, label='2 × monolayer')
ax.set_xlabel(r'Perpendicular bias $F_z$ (V/nm)')
ax.set_xticks(gates)
ax.legend(loc='best',frameon=False,fontsize=6.5)
fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',
           bbox_to_anchor=(.5,.99),ncol=3,frameon=False)
fig.savefig(ROOT / 'MoS2_NMEE_comparison.svg')
svg = ROOT / 'MoS2_NMEE_comparison.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
fig.savefig(ROOT / 'MoS2_NMEE_comparison.pdf')
fig.savefig(ROOT / 'MoS2_NMEE_comparison.png', dpi=600)
plt.close(fig)
