from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
rows=list(csv.DictReader((ROOT/'layer_changes.csv').open()))
fig,axes=plt.subplots(1,2,figsize=(7.2047244094,3.0708661417),layout='constrained')
colors=['#377DA4','#BD6D40'];x=np.arange(2)
for layer,offset,color in zip(['lower','upper'],[-.18,.18],colors):
 r=[next(q for q in rows if q['material']==m and q['layer']==layer) for m in ['MoSSe','InSe']]
 for ax,key,sign in [(axes[0],'G',1),(axes[1],'amplitude_change_percent',-1)]:
  vals=[sign*float(q[key]) for q in r]
  bars=ax.bar(x+offset,vals,.33,color=color,label=layer.capitalize())
  ax.bar_label(bars,labels=[f'{v:.2f}' for v in vals],padding=3,fontsize=8)
for ax,tag in zip(axes,['a','b']):
 ax.set_xticks(x,['MoSSe S/S','InSe A']);ax.axhline(0,color='.4',lw=.6)
 ax.text(-.13,1.05,tag,transform=ax.transAxes,fontweight='bold',fontsize=10)
axes[0].set_ylabel(r'Signed gate slope [$10^{-9}$ (V$^2$ T)$^{-1}$]');axes[0].set_ylim(-14,10)
axes[0].set_title('Absolute layer change',fontsize=9,pad=10)
axes[1].set_ylabel('Decrease of layer amplitude (%)');axes[1].set_ylim(0,5.7)
axes[1].set_title('Relative to each layer’s own response',fontsize=9,pad=10)
axes[1].legend(loc='upper right')
fig.savefig(ROOT/'layer_gate_comparison.svg')
fig.savefig(ROOT/'layer_gate_comparison.pdf')
fig.savefig(ROOT/'layer_gate_comparison.png',dpi=600)
plt.close(fig)

# Normalize generated SVG whitespace without changing geometry.
svg=ROOT/'layer_gate_comparison.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
