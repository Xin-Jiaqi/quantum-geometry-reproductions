"""Quantitative figure from completed independent SCFs only.

Contract: test whether a positive gate changes opposite local-spin magnitudes
in opposite ways in mirror domains, while their net changes have the same sign.
Panel a is the primary local-response evidence; b shows why net tuning alone
misses this distinction. Quantitative grid, 183 x 86 mm, editable PDF/SVG and
600-dpi PNG. No fits, statistical errors, or generated observations.
Highest cutoff/standard probe is displayed; all cutoffs and half-probe checks
remain in source_data.csv, positive_gate_changes.csv and summary.json.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'summary.json').read_text())
rows=s['results']; changes=s['positive_gate_changes']
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'svg.fonttype':'none','pdf.fonttype':42,'savefig.dpi':600})
fig,axs=plt.subplots(1,2,figsize=(183/25.4,86/25.4),gridspec_kw={'width_ratios':[1.12,1]})
fig.subplots_adjust(left=.105,right=.98,bottom=.23,top=.80,wspace=.42)
a,b=axs
for ax,label in zip(axs,'ab'):
    ax.text(-.16,1.19,label,transform=ax.transAxes,fontweight='bold',fontsize=10)
    ax.axhline(0,color='#777777',lw=.7,zorder=0)
colors=['#346995','#BC702C']
for j,layer in enumerate(('lower','upper')):
    for di,domain in enumerate(('A','B')):
        x=di+(j-.5)*.30
        rr=[d for d in changes if d['domain']==domain and d['cutoff']==160 and d['step']=='standard']
        for mesh in (8,6):
            y=next(d[f'delta_abs_C_{layer}'] for d in rr if d['mesh']==mesh)
            a.scatter(x+(-.04 if mesh==6 else .04),y,s=34 if mesh==8 else 25,
                      marker='o' if mesh==8 else 's',facecolor=colors[j] if mesh==8 else 'white',
                      edgecolor=colors[j],linewidth=.9,zorder=4)
a.set_xticks([0,1],['Stack A','Mirror stack B']);a.set_xlim(-.45,1.45)
a.set_ylabel(r'$|C_\ell(F)|-|C_\ell(0)|$'+'\n'+r'[$10^{-9}$ (nm V T)$^{-1}$]')
a.set_title('Local spin-response magnitude',pad=14)
a.legend(handles=[Line2D([],[],color=c,marker='o',lw=0,label=l) for c,l in zip(colors,['Lower layer','Upper layer'])],
         loc='upper left',ncol=1,handletextpad=.4,frameon=False)
a.text(.5,-.28,r'Same positive gate: $F=+0.10$ V/nm',ha='center',transform=a.transAxes,fontsize=8)
for domain,marker,ls,col in [('A','o','-','#343434'),('B','s','--','#777777')]:
    for mesh in (6,8):
        points=sorted([r for r in rows if r['domain']==domain and r['mesh']==mesh and r['cutoff']==160 and r['step']=='standard'],key=lambda r:r['F_V_nm'])
        zero=next(r['C_net'] for r in points if r['F_V_nm']==0)
        xx=[p['F_V_nm'] for p in points];yy=[p['C_net']-zero for p in points]
        if mesh==6:
            b.plot(xx,yy,color=col,ls=ls,lw=1,marker=marker,mfc='white',mec=col,ms=5,label=f'Stack {domain}')
        else:
            b.scatter(xx,yy,marker=marker,color=col,s=23,zorder=5)
b.set_xlabel(r'Perpendicular gate $F$ (V/nm)')
b.set_ylabel(r'$C_{\mathrm{net}}(F)-C_{\mathrm{net}}(0)$'+'\n'+r'[$10^{-9}$ (nm V T)$^{-1}$]')
b.set_xticks([-.1,0,.1]);b.set_xlim(-.115,.115)
b.set_title('Net spin-response change',pad=14)
b.legend(loc='upper right',frameon=False)
fig.text(.54,.96,r'InSe bilayer: $s_{x,\ell}=C_\ell E_x B_z$',ha='center',fontsize=10)
fig.text(.54,.02,'Filled: denser mesh     Open: initial mesh     Fixed ions; lines guide the eye',ha='center',fontsize=7,color='#555555')
for ext in ('svg','pdf','png'):
    fig.savefig(ROOT/f'inse_mirror_gate.{ext}',dpi=600)
plt.close(fig)
print('Saved inse_mirror_gate.svg/.pdf/.png')
