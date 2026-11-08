"""Material comparison: the same lower layer with and without the interface.
Contract: all three independently relaxed strain points and all retained-setting
variants are included. Filled points select an actual calculation, vertical
ranges show setting sensitivity, not statistical error. Straight lines guide
the eye, without interpolation claims. Python, editable SVG/PDF, 89 x 78 mm.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
with (root/'responses.csv').open() as f:rows=list(csv.DictReader(f))
assert set(float(r['strain_pct']) for r in rows)=={-2.,0.,2.}
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'], 'font.size':7.3,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
fig,ax=plt.subplots(figsize=(3.50394,3.07087));fig.subplots_adjust(left=.19,right=.97,bottom=.18,top=.76)
styles=[('A','Bilayer: lower layer','#703AA0','o'),('mono','Matched isolated layer','#4E6372','s'),('electrostatic','Electrostatic reference','#168579','D')]
for kind,label,color,marker in styles:
 x=np.array([-2.,0.,2.]);y=[];lo=[];hi=[]
 for strain in x:
  rr=[r for r in rows if r['system']==kind and float(r['strain_pct'])==strain]
  rep=next(r for r in rr if int(r['nb'])==(96 if kind=='A' else 64) and (kind!='electrostatic' or float(r['reference_end_A'])==3.2))
  yy=float(rep['lower']);y.append(yy);values=[float(r['lower']) for r in rr];lo.append(min(values));hi.append(max(values))
 ax.vlines(x,lo,hi,color=color,lw=2.6,alpha=.35,zorder=2)
 ax.plot(x,y,color=color,marker=marker,ms=5,lw=1.35,label=label,zorder=3)
ax.axhline(0,color='.35',lw=.8,ls='--',zorder=0)
ax.set_xticks([-2,0,2],['−2','0','+2']);ax.set_xlim(-2.25,2.25)
ax.set_xlabel('Biaxial strain (%)')
ax.set_ylabel(r'$C_{xx}^{L}$  [$10^{-9}$ (nm V T)$^{-1}$]')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),borderaxespad=0,fontsize=7.0,handlelength=2.4,labelspacing=.35)
ax.tick_params(direction='out',length=3,width=.8)
fig.savefig(root/'MoSSe_strain_layer_response.svg')
fig.savefig(root/'MoSSe_strain_layer_response.pdf')
fig.savefig(root/'MoSSe_strain_layer_response.png',dpi=600)
plt.close(fig)
