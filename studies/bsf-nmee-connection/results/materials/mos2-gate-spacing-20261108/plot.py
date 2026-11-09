"""Actual fixed-registry reference points only; no fitted sign boundary."""
from pathlib import Path
import csv,json
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
root=Path(__file__).resolve().parent
jobs=json.loads((root/'jobs.json').read_text())
data=[json.loads((root/'baseline.json').read_text())]+[json.loads((root/c/'result.json').read_text()) for c in jobs]
fig,axes=plt.subplots(1,3,figsize=(7.00787,2.75591))
fig.subplots_adjust(left=.09,right=.99,bottom=.29,top=.78,wspace=.35)
fields=sorted({r['physical_Fz_Vnm'] for r in data});colors={0.:'#555555',-.5:'#306A9D',-.25:'#BC793F'}
sources=[]
for field in fields:
 group=sorted([r for r in data if r['physical_Fz_Vnm']==field],key=lambda r:r['delta_d_A'])
 x=[r['delta_d_A'] for r in group];color=colors.get(field,'#775B87')
 for nb in ('112','96'):
  for ax,comp in zip(axes,(1,2,0)):
   vals=[r['results'][nb]['standard']['C_total_lower_upper'][comp] for r in group]
   ax.plot(x,vals,marker='s' if field==0 else 'o',color=color,lw=1.1 if nb=='96' else .6,
           alpha=1 if nb=='96' else .35,ms=4,mfc=color if nb=='96' else 'white',
           label=f'{field:+.2f} V/nm' if nb=='96' else None)
   for r,v in zip(group,vals):sources.append({'case':r['case'],'delta_d_A':r['delta_d_A'],'Fz_Vnm':field,
       'states':nb,'component':{0:'total',1:'lower',2:'upper'}[comp],'C':v})
for ax,title,letter in zip(axes,('Lower layer','Upper layer','Net response'),('a','b','c')):
 ax.axhline(0,color='.45',lw=.7,ls=':');ax.set_title(title,fontsize=9)
 ax.set_xlabel(r'Separation increase $\Delta d$ ($\mathrm{\AA}$)')
 ax.set_xticks([0,.15,.3]);ax.text(-.18,1.13,letter,weight='bold',transform=ax.transAxes)
axes[0].set_ylabel(r'$C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=len(fields),bbox_to_anchor=(.52,1.015))
fig.text(.5,.055,'Fixed 3R registry; faint curves show retained-state sensitivity. Lines guide the eye.',ha='center',fontsize=7,color='.35')
fig.savefig(root/'MoS2_gate_spacing.png',dpi=600)
fig.savefig(root/'MoS2_gate_spacing.pdf')
fig.savefig(root/'MoS2_gate_spacing.svg')
svg=root/'MoS2_gate_spacing.svg';svg.write_text('\n'.join(x.rstrip() for x in svg.read_text().splitlines())+'\n')
with (root/'figure_source.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(sources[0]),lineterminator='\n');w.writeheader();w.writerows(sources)
print('plotted',len(data),'actual points;',len(sources),'component/cutoff values')
