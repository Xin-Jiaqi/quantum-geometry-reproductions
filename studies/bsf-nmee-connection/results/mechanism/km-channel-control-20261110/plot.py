from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
p=json.loads((ROOT/'weak_prediction.json').read_text());s=json.loads((ROOT/'summary.json').read_text())
fig,(ax,bx)=plt.subplots(1,2,figsize=(7.2047244094,3.1496062992),layout='constrained')
alpha=np.linspace(-1,1,201);a,b,c=[p[k][1]*1e3 for k in ['A_background','B_odd_channel_interference','D_even_channel']]
ax.plot(alpha,a+b*alpha+c*alpha**2,color='#BD6D40',lw=1.5,label='Weak-coupling prediction')
rr=json.loads((ROOT/'weak_training.json').read_text())['rows'];ax.scatter([r['retained_alpha'] for r in rr],[r['G'][1]*1e3 for r in rr],s=23,c='#BD6D40',marker='s',label='Coefficient inputs',zorder=4)
vv=[r for r in s['rows'] if r['t_meV']==2 and r['retained_alpha'] in [-.5,.5] and r['source'] in ['weak_validation.json','phase_validation.json']]
ax.scatter([r['retained_alpha'] for r in vv],[r['G_upper']*1e3 for r in vv],s=48,facecolors='white',edgecolors='#222222',linewidths=1.2,label='Held-out calculations',zorder=5)
ax.axhline(0,color='.5',lw=.7);ax.set_xlabel(r'Signed selected-channel amplitude $\alpha$');ax.set_ylabel(r'$G_\mathrm{upper}$ [$10^{-3}$ (nm V T eV)$^{-1}$]');ax.set_title('Gate direction depends on interference',fontsize=9)
ax.legend(loc='lower center',fontsize=7);ax.set_ylim(-3.2,.95)
cs=s['strong_contrast'];x=np.arange(2)
for offset,key,color,label in [(-.18,'G_upper_channel_removed','#8194A1','Selected channel removed'),(.18,'G_upper_full','#BD6D40','Full interface')]:
 vals=[abs(r[key])*1e3 for r in cs];bars=bx.bar(x+offset,vals,.33,color=color,label=label);bx.bar_label(bars,labels=[f'{v:.1f}' for v in vals],padding=3,fontsize=8)
bx.set_xticks(x,['10 meV','20 meV']);bx.set_xlabel(r'Overall interlayer hopping $t_\perp$');bx.set_ylabel(r'$|G_\mathrm{upper}|$ [$10^{-3}$ (nm V T eV)$^{-1}$]');bx.set_ylim(0,80);bx.legend(loc='upper left',fontsize=7);bx.set_title('The same channel can enhance or suppress',fontsize=9)
for a1,letter in [(ax,'a'),(bx,'b')]:a1.text(-.13,1.05,letter,transform=a1.transAxes,fontweight='bold',fontsize=10)
fig.savefig(ROOT/'km_channel_mechanism.svg');fig.savefig(ROOT/'km_channel_mechanism.pdf');fig.savefig(ROOT/'km_channel_mechanism.png',dpi=600);plt.close(fig)
f=ROOT/'km_channel_mechanism.svg';f.write_text('\n'.join(x.rstrip() for x in f.read_text().splitlines())+'\n')
