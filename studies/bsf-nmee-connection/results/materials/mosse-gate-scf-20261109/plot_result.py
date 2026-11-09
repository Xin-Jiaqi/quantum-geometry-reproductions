"""Screening test: absolute net gate sensitivity and signed layer sensitivities.

Quantitative two-panel comparison. Panel a tests control gain; panel b tests
upper-layer dominance independently. Every tested state cutoff is retained.
No statistical error interpretation or near-zero-baseline enhancement ratio.
The same ±0.10 V/nm interval is used for model and self-consistent derivatives.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
s=json.loads((root/'summary.json').read_text());assert s['complete']
m=pd.read_csv(root/'model_slopes.csv');m=m[np.isclose(m.half_interval,.1)]
m=m[(m.kind=='isolated_pair')|((m.kind=='full')&(m.frozen_empty==8))]
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':600})
fig,axs=plt.subplots(1,2,figsize=(7.2047244,2.9527559),layout='constrained')
sources=[]
def series(ax,x,vals,color,marker,label=None):
    vals=np.array(vals)
    ax.plot([x,x],[vals.min(),vals.max()],color=color,lw=1.4,zorder=2)
    ax.scatter(np.full(len(vals),x),vals,s=24,marker=marker,facecolor='white',
               edgecolor=color,lw=.9,label=label,zorder=3)
    ax.scatter([x],[vals[-1]],s=25,marker=marker,color=color,zorder=4)

for offset,kind,color,marker,label in [(-.13,'isolated_pair','#6B7B8A','s','Isolated-layer pair'),
                                     (.13,'full','#7344A4','o','Bilayer')]:
    vals0=abs(m.loc[m.kind==kind,'G_net'].values)
    if kind=='full':vals1=[abs(s['slopes']['A'][nb][0]) for nb in ('96','112')]
    else:vals1=[abs(s['slopes']['mono'][nb][0]+s['slopes']['upper'][nb][0]) for nb in ('48','56','64')]
    for x,vals,condition in [(0,vals0,'model'),(1,vals1,'SCF')]:
        series(axs[0],x+offset,vals,color,marker,label if x==0 else None)
        sources.extend({'panel':'a','condition':condition,'kind':kind,'value':float(v)} for v in vals)
axs[0].set_xticks([0,1],['Gate screening\nomitted','Gate screening\nincluded'])
axs[0].set_ylabel(r'$|\Delta C_{\rm net}/\Delta F|$  [$10^{-9}$ (V$^2$ T)$^{-1}$]')
axs[0].set_ylim(bottom=0);axs[0].set_xlim(-.45,1.45)
axs[0].legend(frameon=False,loc='upper right')
for off,component,color,marker,label in [(-.15,1,'#356FB2','s','Lower'),
                                       (0,2,'#C75B38','^','Upper'),(.15,0,'#252525','o','Net')]:
    iso=[]
    for nb in ('48','56','64'):
        gl=s['slopes']['mono'][nb][0];gu=s['slopes']['upper'][nb][0]
        iso.append([gl+gu,gl,gu][component])
    full=[s['slopes']['A'][nb][component] for nb in ('96','112')]
    for x,vals,kind in [(0,iso,'isolated_pair'),(1,full,'full')]:
        series(axs[1],x+off,vals,color,marker,label if x==0 else None)
        sources.extend({'panel':'b','condition':'SCF','kind':kind+'_'+label,'value':float(v)} for v in vals)
axs[1].set_xticks([0,1],['Isolated layers','Bilayer'])
axs[1].set_ylabel(r'$\Delta C/\Delta F$  [$10^{-9}$ (V$^2$ T)$^{-1}$]')
axs[1].axhline(0,color='.7',lw=.6);axs[1].set_xlim(-.5,1.5)
axs[1].margins(y=.23);axs[1].legend(frameon=False,ncol=3,loc='upper center',
 handletextpad=.3,columnspacing=.7)
for i,ax in enumerate(axs):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3,width=.6)
    ax.text(-.14,1.02,'ab'[i],transform=ax.transAxes,fontweight='bold',fontsize=10)
fig.savefig(root/'gate_scf.png',dpi=600)
fig.savefig(root/'gate_scf.pdf')
fig.savefig(root/'gate_scf.svg')
p=root/'gate_scf.svg';p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
pd.DataFrame(sources).to_csv(root/'figure_source.csv',index=False)
(root/'figure_contract.json').write_text(json.dumps({'question':'Does gate screening preserve net gain and upper-layer dominance?',
 'archetype':'quantitative grid','backend':'python','dimensions_mm':[183,75],
 'data':'same ±0.10 interval; all tested state cutoffs retained, duplicate full projection windows omitted',
 'mark_definition':'points and spans: basis choices; filled point highest state cutoff, not statistical uncertainty',
 'panel_a':'absolute net sensitivity, independent layer sum vs bilayer, with/without gate screening',
 'panel_b':'signed layer and net sensitivities after gate SCF',
 'no_near_zero_normalization':True,'model_upper_reference':'approximate mirror; SCF upper reference exact geometry'},indent=2))
