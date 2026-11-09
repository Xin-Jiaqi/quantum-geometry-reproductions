from pathlib import Path
import json,csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
study=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent/'figures'
out.mkdir(exist_ok=True)
r=study/'results/mechanism'
f=json.load(open(study/'docs/revision-20261106c/figure_data.json'))
sets={n:json.load(open(r/'zero-field-20261107/refined'/name)) for n,name in [(96,'controls96-analysis.json'),(112,'full-controls-analysis.json')]}
channels=json.load(open(r/'mixing-channels-20261107/material-channels.json'))
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
blue='#315B88';orange='#BE793C';purple='#7030A0';gray='#B5C0CC'
rows=[]
def save(fig,name):
 fig.savefig(out/(name+'.png'),dpi=600)
 fig.savefig(out/(name+'.pdf'))
 fig.savefig(out/(name+'.svg'))
 svg=out/(name+'.svg')
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
 plt.close(fig)
def zero(ax):ax.axhline(0,color='.4',lw=.65,ls=':')
fig,ax=plt.subplots(1,3,figsize=(7.1,2.65));fig.subplots_adjust(left=.085,right=.985,bottom=.25,top=.84,wspace=.43)
for n,d in sets.items():
 for w in ('4','8','12'):
  rr=d['windows'][w]['hopping_scan'];main=n==96 and w=='8'
  ax[0].plot([x['alpha'] for x in rr],[x['C'] for x in rr],'o-',color=blue if main else gray,lw=1.5 if main else .65,ms=3 if main else 2,zorder=4 if main else 1)
  rows += [dict(figure='7',panel='a',series=f'{n}_window{w}',x=x['alpha'],y=x['C']) for x in rr]
ax[0].set(xlabel=r'Interlayer mixing fraction $\alpha$',ylabel=r'$C_{xx}$  [$10^{-9}$ (nm V T)$^{-1}$]')
ax[0].set_title('Same lower layer',fontsize=9)
pairs=f['mirror']['paired'];xs=[p['F_A'] for p in pairs]
for label,idx,c in [('Lower',1,blue),('Upper',2,orange),('Total',0,'#222222')]:
 ys=[p['A'][idx] for p in pairs];ax[1].plot(xs,ys,'o-',color=c,lw=1,ms=3,label=label)
 rows += [dict(figure='7',panel='b',series=label,x=x,y=y) for x,y in zip(xs,ys)]
ax[1].set(xlabel=r'Control field $F_z$ (V/nm)');ax[1].legend(fontsize=7,loc='upper left',bbox_to_anchor=(.02,.78));ax[1].set_title('Fixed S/S stacking A',fontsize=9)
for tag,c,sgn in [('A',blue,1),('B',orange,-1)]:
 xx=[sgn*p['F_A'] for p in pairs];yy=[p[tag][0] for p in pairs];sort=np.argsort(xx)
 ax[2].plot(np.array(xx)[sort],np.array(yy)[sort],'o-',c=c,lw=1,ms=3,label=tag)
 rows += [dict(figure='7',panel='c',series=tag,x=x,y=y) for x,y in zip(xx,yy)]
ax[2].set(xlabel=r'Control field $F_z$ (V/nm)');ax[2].legend(fontsize=7);ax[2].set_title('Mirror-related stackings',fontsize=9)
for i,a in enumerate(ax):zero(a);a.text(-.20,1.13,'abc'[i],weight='bold',transform=a.transAxes)
fig.text(.5,.025,'(a) Zero-field material-model control.  (b,c) Independent finite-field DFT; original field convention.',ha='center',fontsize=7,color='.35')
save(fig,'MoSSe_controls_zero_reference')
fig,ax=plt.subplots(1,2,figsize=(7.1,3.05),gridspec_kw={'width_ratios':[1,1.13]});fig.subplots_adjust(left=.085,right=.985,bottom=.31,top=.84,wspace=.27)
keys=['states0_energy0','states0_energy1','states1_energy0','states1_energy1']
labels=['Weak states\nWeak energies','Weak states\nFull energies','Full states\nWeak energies','Full states\nFull energies']
for j,key in enumerate(keys):
 vals=[]
 for n,d in sets.items():
  for w in ('4','8','12'):
   v=d['windows'][w]['spectral'][key]['C'];vals.append(v);rows.append(dict(figure='8',panel='a',series=f'{n}_window{w}_{key}',x=j,y=v))
 ax[0].scatter(j+np.linspace(-.07,.07,len(vals)),vals,s=16,facecolors='white',edgecolors=orange if j<2 else blue,lw=.8)
ax[0].set(xticks=range(4),xticklabels=labels,ylabel=r'$C_{xx}^{\rm lower}$  [$10^{-9}$ (nm V T)$^{-1}$]',ylim=(-17,24))
ax[0].set_title('The sign follows the electronic states',fontsize=9)
for j,key in enumerate(['none','VV','CC','VV_CC','VC','all']):
 vals=[]
 for name,d in channels['runs'].items():
  v=d['rows'][key]['response']['C'];vals.append(v);rows.append(dict(figure='8',panel='b',series=name+'_'+key,x=j,y=v))
 ax[1].scatter(j+np.linspace(-.08,.08,len(vals)),vals,s=16,facecolors='white',edgecolors=blue if key in ('VC','all') else orange,lw=.8)
ax[1].set(xticks=range(6),xticklabels=['None','V–V','C–C','V–V\n+ C–C','V–C','Full'],ylim=(-17,24))
ax[1].set_title('Cross-gap mixing is sufficient to reverse',fontsize=9)
for i,a in enumerate(ax):zero(a);a.tick_params(axis='x',labelsize=7,length=0,pad=8);a.text(-.14,1.13,'ab'[i],weight='bold',transform=a.transAxes)
fig.text(.5,.04,'MoSSe S/S, fixed structure and interface potential.  V: occupied; C: empty states.\nEach dot is a material-model setting; channel interventions are not additive contributions.',ha='center',fontsize=7,color='.35',linespacing=1.4)
save(fig,'MoSSe_states_and_channels')
with open(out/'figure_source.csv','w') as stream:
 w=csv.DictWriter(stream,fieldnames=['figure','panel','series','x','y'],lineterminator='\n');w.writeheader();w.writerows(rows)
