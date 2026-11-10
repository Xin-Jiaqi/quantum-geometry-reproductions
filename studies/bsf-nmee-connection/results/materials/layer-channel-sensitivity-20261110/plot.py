from pathlib import Path
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
rows=list(csv.DictReader((ROOT/'all_controls.csv').open()))
fig,axes=plt.subplots(1,2,figsize=(7.2047244094,3.3858267717),layout='constrained',sharey=True)
for ax,mat,letter in zip(axes,['MoSSe','InSe'],['a','b']):
 for i,ch in enumerate(['full','Vb_Ct']):
  sub=[r for r in rows if r['material']==mat and r['channel']==ch]
  representative=next(r for r in sub if r['representative']=='True')
  for layer,offset,color in [('lower',-.18,'#377DA4'),('upper',.18,'#BD6D40')]:
   key='amplitude_change_percent_'+layer;v=float(representative[key]);x=i+offset
   ax.bar(x,v,.33,color=color,label=layer.capitalize() if i==0 else None,zorder=2)
   ax.scatter([x]*len(sub),[float(r[key]) for r in sub],s=9,facecolors='none',edgecolors='.2',linewidths=.5,zorder=3)
   label_x=x-.22 if mat=='MoSSe' and ch=='Vb_Ct' and layer=='lower' else x
   ax.text(label_x,v+(0.22 if v>=0 else -.3),f'{v:+.2f}',ha='center',va='bottom' if v>=0 else 'top',fontsize=8)
 ax.axhline(0,color='.2',lw=.7);ax.set_xticks([0,1],['Full interface',r'Remove $V_\mathrm{lower}$–$C_\mathrm{upper}$'])
 ax.set_title('MoSSe S/S' if mat=='MoSSe' else 'InSe A',fontsize=10)
 ax.text(-.12,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
 ax.set_ylim(-9.3,2.4)
axes[0].set_ylabel('Change of layer response amplitude (%)')
axes[1].legend(loc='lower right')
for ax in axes:ax.tick_params(axis='x',labelsize=8)
fig.savefig(ROOT/'channel_gate_comparison.svg')
fig.savefig(ROOT/'channel_gate_comparison.pdf')
fig.savefig(ROOT/'channel_gate_comparison.png',dpi=600)
plt.close(fig)
p=ROOT/'channel_gate_comparison.svg';p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
