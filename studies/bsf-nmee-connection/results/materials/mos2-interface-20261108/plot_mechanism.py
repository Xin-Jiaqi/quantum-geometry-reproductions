"""Fixed-geometry physical controls, with all tested window variations shown."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
p=Path(__file__).resolve().parent;a=json.loads((p/'analysis.json').read_text());e=json.loads((p/'electrostatic_audit.json').read_text());mix=a['materials']['A']['mixing']
fig,axes=plt.subplots(1,2,figsize=(7.00787,2.91339),gridspec_kw={'width_ratios':[1,1.23]});fig.subplots_adjust(left=.1,right=.98,bottom=.30,top=.86,wspace=.32)
blue='#315B88';orange='#BE793C';gray='#BEC8D2';rows=[]
ax=axes[0]
for key,val in mix.items():
 x=[r['alpha'] for r in val['rows']];y=[r['C_total_lower_upper'][1] for r in val['rows']]
 ax.plot(x,y,'o-',color=blue if key=='96_8' else gray,linewidth=1.5 if key=='96_8' else .7,markersize=4 if key=='96_8' else 2.5,zorder=3 if key=='96_8' else 1)
 for xx,yy in zip(x,y):rows.append({'panel':'a','control':key,'parameter':xx,'C_lower':yy})
ax.axhline(0,color='.35',lw=.65,ls=':');ax.set(xlabel='Interlayer mixing fraction',ylabel=r'$C_{xx}^{\mathrm{lower}}$  [$10^{-9}$ (nm V T)$^{-1}$]',xticks=[0,.5,1],ylim=(-6,5))
ax.set_title('Same layer, opposite response',fontsize=9,pad=9);ax.text(-.19,1.13,'a',weight='bold',transform=ax.transAxes)
ax.text(.04,.08,'Fixed atomic structure',transform=ax.transAxes,fontsize=7,color='.3')
sets=[('Isolated\nmonolayer',[r['C_total_lower_upper'][0] for r in a['materials']['mono']['cutoffs'].values()]),
 ('Interface\nelectrostatics',[v['controls'][f'lower_{end}_full']['C_total_lower_upper'][0] for v in e.values() for end in (3.2,3.8)]),
 ('Selected mixing\nsuppressed',[v['rows'][0]['C_total_lower_upper'][1] for v in mix.values()]),
 ('Full\ninterface',[v['rows'][-1]['C_total_lower_upper'][1] for v in mix.values()])]
ax=axes[1]
for i,(label,vals) in enumerate(sets):
 color=orange if i==1 else blue
 jitter=np.linspace(-.08,.08,len(vals));ax.scatter(i+jitter,vals,s=15,facecolors='white',edgecolors=color,linewidths=.8,zorder=3)
 for n,yy in enumerate(vals):rows.append({'panel':'b','control':label.replace('\n',' '),'parameter':n,'C_lower':yy})
ax.axhline(0,color='.35',lw=.65,ls=':');ax.set(xticks=range(len(sets)),xticklabels=[s[0] for s in sets],ylim=(-6,5),xlim=(-.45,3.45))
ax.tick_params(axis='x',labelsize=7,length=0,pad=8);ax.set_title('Smooth electrostatics: opposite sign',fontsize=9,pad=9)
ax.text(-.12,1.13,'b',weight='bold',transform=ax.transAxes)
fig.text(.5,.025,'Frozen-potential controls; dots show window / embedding sensitivity, not statistical uncertainty.',ha='center',fontsize=7,color='.35')
fig.savefig(p/'MoS2_interface_mechanism.png',dpi=350)
fig.savefig(p/'MoS2_interface_mechanism.svg')
svg=p/'MoS2_interface_mechanism.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
fig.savefig(p/'MoS2_interface_mechanism.pdf')
with (p/'figure_source.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['panel','control','parameter','C_lower'],lineterminator='\n');w.writeheader();w.writerows(rows)
print('Figure written, source rows:',len(rows))
