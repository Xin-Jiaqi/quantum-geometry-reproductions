"""Compare compensation shifts using all actual net-response points.

Quantitative comparison, 178 x 80 mm, Python only. The two panels isolate the
full material from its smooth screened electrostatic representation. Curves
connect actual points; none is an independently calculated phase boundary.
All retained-state/spatial-taper settings are shown, not statistical replicates.
"""
from pathlib import Path
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8, 'axes.spines.top':False, 'axes.spines.right':False,
    'axes.linewidth':.7, 'svg.fonttype':'none', 'pdf.fonttype':42,
    'legend.frameon':False})
full = list(csv.DictReader((ROOT / 'full_response.csv').open()))
static = list(csv.DictReader((ROOT / 'electrostatic_response.csv').open()))
assert all(float(r['C_lower']) < 0 < float(r['C_upper']) for r in full)
assert all(float(r['C_upper']) < 0 < float(r['C_lower']) for r in static)
fig, axes = plt.subplots(1, 2, figsize=(7.00787, 3.14961))
fig.subplots_adjust(left=.10, right=.98, bottom=.27, top=.79, wspace=.35)
colors = {0:'#555555', -.1:'#3274A1', -.2:'#B36535'}
markers = {0:'o', -.1:'s', -.2:'^'}
for distance in (0, -.1, -.2):
    for cutoff in (96, 112):
        rows = sorted([r for r in full if float(r['delta_d_A']) == distance and int(r['cutoff']) == cutoff],
                      key=lambda r:float(r['Fz_Vnm']))
        axes[0].plot([float(r['Fz_Vnm']) for r in rows], [float(r['C_total']) for r in rows],
            color=colors[distance], marker=markers[distance], ms=4,
            lw=1.15 if cutoff == 96 else .7, ls='-' if cutoff == 96 else '--',
            alpha=1 if cutoff == 96 else .55, mfc=colors[distance] if cutoff == 96 else 'white',
            label=f'{-distance:.2f} '+r'$\mathrm{\AA}$' if cutoff == 96 else None)
    for cutoff in (40, 48, 56, 64):
        for taper in (3.2, 3.8):
            rows = sorted([r for r in static if float(r['delta_d_A']) == distance and int(r['cutoff']) == cutoff
                           and float(r['taper_end_A']) == taper], key=lambda r:float(r['Fz_Vnm']))
            axes[1].plot([float(r['Fz_Vnm']) for r in rows], [float(r['C_total']) for r in rows],
                color=colors[distance], marker=markers[distance], ms=2.7, lw=.7,
                alpha=.42, mfc='white', ls='-' if taper == 3.2 else '--')
for ax, label, title in zip(axes, ('a', 'b'), ('Full bilayer', 'Screened electrostatic reference')):
    ax.axhline(0, color='.35', lw=.75, ls=':')
    ax.set_xlim(-.52, .02); ax.set_xticks([-.5, -.32, 0])
    ax.set_xlabel(r'Control field $F_z$ (V nm$^{-1}$)')
    ax.set_ylabel(r'Net $C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]')
    ax.set_title(title, fontsize=9)
    ax.text(-.15, 1.11, label, transform=ax.transAxes, weight='bold', fontsize=10)
axes[0].text(.04, .10, 'Layer directions: lower −, upper +', transform=axes[0].transAxes, fontsize=7.5)
axes[1].text(.04, .10, 'Layer directions: lower +, upper −', transform=axes[1].transAxes, fontsize=7.5)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, title='Compression of the interlayer gap', loc='upper center',
           bbox_to_anchor=(.53, 1.02), ncol=3, fontsize=8, title_fontsize=8)
fig.text(.5, .07, 'Actual field points; lines guide the eye. Curve spread is model sensitivity, not statistical error.',
         ha='center', fontsize=7, color='.3')
fig.savefig(ROOT / 'MoS2_compression_compensation.png', dpi=600)
fig.savefig(ROOT / 'MoS2_compression_compensation.pdf')
fig.savefig(ROOT / 'MoS2_compression_compensation.svg')
svg = ROOT / 'MoS2_compression_compensation.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
print('Plotted all', len(full), 'full and', len(static), 'electrostatic net-response values.')
