"""Plot calculated gate increments and layer sensitivity, without fitting.

Contract: interface mixing moderately enhances gate sensitivity and makes the
upper layer dominate. Panel a compares absolute net increments, avoiding a
near-zero denominator. Panel b shows signed layer slopes. Thin lines/ranges are
model-setting variation, not statistical uncertainty. Full curves duplicated
across projection windows are deduplicated because no projection is applied.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
df=pd.read_csv(root/'curves.csv');sl=pd.read_csv(root/'slopes.csv')
sl=sl[np.isclose(sl.half_interval,.025)]
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
    'axes.titlesize':8,'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7,
    'axes.linewidth':.65,'lines.linewidth':1.3,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','savefig.dpi':600})
fig,axs=plt.subplots(1,2,figsize=(7.2047244,2.8740157),layout='constrained',
                     gridspec_kw={'width_ratios':[1.05,1]})
colors={'full':'#6B40A0','electrostatic_pair':'#D48627','isolated_pair':'#3979A3',
        'no_VC':'#5E7167'}
labels={'full':'Full bilayer','electrostatic_pair':'Electrostatic reference',
        'isolated_pair':'Isolated-layer pair','no_VC':'Cross-layer V–C removed'}
ax=axs[0]
for kind in ('full','electrostatic_pair','isolated_pair'):
    dat=df[df.kind==kind]
    if kind=='full':dat=dat[dat.frozen_empty==8]
    for _,g in dat.groupby(['nb','frozen_empty','taper_end_A'],dropna=False):
        g=g.sort_values('F_model');baseline=g.loc[np.isclose(g.F_model,0),'C_net'].iloc[0]
        ax.plot(g.F_model,g.C_net-baseline,color=colors[kind],alpha=.23,lw=.75)
    if kind=='full':main=dat[dat.nb==112]
    elif kind=='electrostatic_pair':main=dat[(dat.nb==64)&np.isclose(dat.taper_end_A,3.2)]
    else:main=dat[dat.nb==64]
    main=main.sort_values('F_model');base=main.loc[np.isclose(main.F_model,0),'C_net'].iloc[0]
    ax.plot(main.F_model,main.C_net-base,'o-',ms=2.5,color=colors[kind],label=labels[kind])
ax.axhline(0,color='.65',lw=.5,zorder=0);ax.axvline(0,color='.65',lw=.5,zorder=0)
ax.set(xlabel=r'Model gate field $F$ (V nm$^{-1}$)',
       ylabel=r'$C_{\rm net}(F)-C_{\rm net}(0)$  [$10^{-9}$ (nm V T)$^{-1}$]')
ax.legend(frameon=False,loc='upper right',handlelength=1.8)
ax.text(-.14,1.02,'a',transform=ax.transAxes,fontweight='bold',fontsize=10)
ax=axs[1];kinds=['isolated_pair','electrostatic_pair','no_VC','full']
for offset,col,marker,color,label in [(-.15,'G_lower','s','#316BB0','Lower'),
        (0,'G_upper','^','#C95E3D','Upper'),(.15,'G_net','o','#252525','Net')]:
    for x,kind in enumerate(kinds):
        vals=sl.loc[sl.kind==kind,col].values
        ax.plot([x+offset]*2,[vals.min(),vals.max()],color=color,lw=1.5)
        ax.scatter([x+offset],[np.mean([vals.min(),vals.max()])],marker=marker,
                   s=20,color=color,label=label if x==0 else None,zorder=3)
ax.set_xticks(range(4),['Isolated\nlayers','Electrostatic\nreference','V–C\nremoved','Full\nbilayer'])
ax.set(ylabel=r'$\partial C/\partial F$  [$10^{-9}$ (V$^2$ T)$^{-1}$]',ylim=(-18,0.7))
ax.axhline(0,color='.65',lw=.5);ax.legend(frameon=False,ncol=3,loc='upper right',
                                     handletextpad=.3,columnspacing=.8)
ax.text(-.17,1.02,'b',transform=ax.transAxes,fontweight='bold',fontsize=10)
for ax in axs:
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(direction='out',length=3,width=.6)
fig.savefig(root/'gate_efficiency.png',dpi=600)
fig.savefig(root/'gate_efficiency.pdf')
fig.savefig(root/'gate_efficiency.svg')
svg=root/'gate_efficiency.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
(root/'figure_contract.json').write_text(json.dumps({
    'conclusion':'Moderate net gate-sensitivity enhancement, dominated by upper-layer response.',
    'archetype':'comparison plus layer-resolved mechanism control',
    'data':['curves.csv','slopes.csv'],'no_fitting':True,
    'panel_a':'all computed increments; duplicate full curves removed; high-band-count representative emphasized',
    'panel_b':'signed central derivatives at delta F=0.025; marker midpoint and span of settings, not error bars',
    'no_statistical_interpretation':True,'dimensions_mm':[183,73],
    'field_scope':'one-electron imposed slope, not calibrated external voltage',
    'exports':['PDF','SVG','PNG 600 dpi']},indent=2))
