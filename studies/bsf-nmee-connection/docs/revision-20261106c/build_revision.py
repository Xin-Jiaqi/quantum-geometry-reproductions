"""Assemble a bounded KM-registry revision from archived response calculations."""
from pathlib import Path
import csv
import json
import numpy as np

R = Path(__file__).resolve().parent
B = R.parent / 'revision-20261106b'


def read_rows(name):
    with (R / name).open() as f:
        rows = list(csv.DictReader(f))
    return [{k: (v if k == 'grid' else float(v)) for k, v in row.items()} for row in rows]


D = json.loads((B / 'figure_data.json').read_text())
D['km_normal'] = read_rows('normal_main.csv') + read_rows('gate_extra_main.csv')
D['km_neighbour'] = read_rows('neighbourhood_robust.csv') + read_rows('gate_extra_robust.csv')
D['km_tensor_audit'] = json.loads((R / 'tensor_check.json').read_text())
(R / 'figure_data.json').write_text(json.dumps(D, indent=2, ensure_ascii=False) + '\n')

# Preserve the previous AB/BA figure as a supplementary control.
script = (B / 'plot_revision.py').read_text()
script = script.replace('figures-20261106b', 'figures-20261106c')
script = script.replace("save(fig,'Fig1_KM')", "save(fig,'FigS4_KM_AB_BA')")
script = script.replace("'AB and its mirror partner'", "'KM AB and its mirror partner'")
script = script.replace("'Hybridization at AB'", "'Hybridization at KM AB'")
script = script.replace("'Layer detuning at AB'", "'Layer detuning at KM AB'")
script = script.replace("xticklabels=['AB','BA']", "xticklabels=['KM AB','KM BA']")
start = script.index('# The previously displayed generic displacement')
end = script.index('\ndef material', start)
script = script[:start] + r'''# New main figure: symmetry-defined probe at a fixed intermediate KM registry.
def quarter_geometry(ax):
 lattice=np.array([[-np.sqrt(3)/2,1.5],[np.sqrt(3)/2,1.5]])
 base=np.array([np.array([i,j])@lattice+q for i in range(-3,4) for j in range(-3,4) for q in ([0,0],[0,1])])
 shift=np.array([-1/6,7/12])@lattice
 # Display x along the mirror normal; y completes a right-handed basis.
 rot=np.array([[.5,np.sqrt(3)/2],[-np.sqrt(3)/2,.5]])
 for off,marker,fc,col in [(np.zeros(2),'o',T,T),(shift,'s','white',B)]:
  xy=(base+off)@rot.T
  mask=(abs(xy[:,0])<2.15)&(abs(xy[:,1])<2.15);xy=xy[mask]
  for i in range(len(xy)):
   for j in range(i):
    if abs(np.linalg.norm(xy[i]-xy[j])-1)<1e-6:ax.plot(xy[[i,j],0],xy[[i,j],1],c=col,lw=.6,alpha=.5)
  ax.scatter(xy[:,0],xy[:,1],s=24,marker=marker,facecolors=fc,edgecolors=col,lw=.9,zorder=3)
 ax.axvline(0,ymin=.20,ymax=.88,c=A,ls='--',lw=.9,zorder=0)
 ax.text(.12,2.13,r'$m_v$',c=A,fontsize=8)
 ax.annotate('',xy=(1.75,-2.5),xytext=(.3,-2.5),arrowprops={'arrowstyle':'->','color':N,'lw':1.2})
 ax.text(1.1,-2.95,r'$E_x$',ha='center',fontsize=8)
 ax.text(-1.8,-2.65,r'$\odot\ B_z$',fontsize=9)
 ax.set(xlim=(-2.4,2.4),ylim=(-3.15,2.6));ax.set_aspect('equal');ax.axis('off')
 ax.text(.02,1.01,'Upper',color=T,transform=ax.transAxes,fontsize=7)
 ax.text(.68,1.01,'Lower',color=B,transform=ax.transAxes,fontsize=7)
 ax.text(.5,-.08,r'KM AB  $\longrightarrow$  KM BA',transform=ax.transAxes,ha='center',fontsize=7)
 ax.text(.5,-.18,r'Fixed registry: $\lambda=1/4$',transform=ax.transAxes,ha='center',fontsize=8)

fig,axs=plt.subplots(2,2,figsize=(7,5.8));fig.subplots_adjust(left=.11,right=.97,bottom=.10,top=.91,wspace=.40,hspace=.68)
a,b,c,d=axs.flat;quarter_geometry(a);panel(a,'a','Mirror-defined probing direction')
r=sorted([v for v in D['km_normal'] if v['U_top_minus_bottom_meV']==0],key=lambda v:v['t_meV'])
layers(b,r,'t_meV',1e-5);panel(b,'b','Reversing the lower-layer spin')
b.set(xlabel=r'$t_\perp$ (meV)',ylabel=unit5);b.legend(loc='upper left',fontsize=7)
b.text(.97,.04,r'$U=0$',transform=b.transAxes,ha='right',fontsize=8)
r2=sorted([v for v in D['km_normal'] if v['t_meV']==10 and v['U_top_minus_bottom_meV']<=5],key=lambda v:v['U_top_minus_bottom_meV'])
layers(c,r2,'U_top_minus_bottom_meV',1e-5);panel(c,'c','Reversing only the total spin')
c.set(xlabel=r'$U$ (meV)',ylabel=unit5,xlim=(-.1,5.1))
c.text(.97,.96,r'$t_\perp=10$ meV',transform=c.transAxes,ha='right',va='top',fontsize=8)
for lam,col,marker in [(.24,A,'o'),(.25,N,'s'),(.26,P,'^')]:
 rr=sorted([v for v in D['km_neighbour'] if v['mass_meV']==45 and v['lambda']==lam and v['t_meV']==10 and 3<=v['U_top_minus_bottom_meV']<=5],key=lambda v:v['U_top_minus_bottom_meV'])
 d.plot([v['U_top_minus_bottom_meV'] for v in rr],[v['total']*1e5 for v in rr],marker+'-',c=col,ms=3,lw=1.1,label=rf'$\lambda={lam}$')
zero(d);panel(d,'d','Nearby fixed KM registries');d.set(xlabel=r'$U$ (meV)',ylabel=r'Total '+unit5)
d.legend(loc='upper right',fontsize=7);d.text(.04,.05,r'$t_\perp=10$ meV',transform=d.transAxes,fontsize=7)
save(fig,'Fig1_KM')

# Disclose the old projection zero and the continuation beyond the new first zero.
fig,(a,b)=plt.subplots(1,2,figsize=(7,3.05));fig.subplots_adjust(left=.10,right=.97,bottom=.20,top=.80,wspace=.40)
rr=sorted([v for v in D['km_tensor_audit'] if v['lambda']==.25],key=lambda v:v['U_meV'])
for i,col,lab in [(0,N,r'$C_{XX}^{\mathrm{tot}}$'),(1,A,r'$C_{YX}^{\mathrm{tot}}$')]:
 a.plot([v['U_meV'] for v in rr],[v['tensor']['total'][i][0]*1e5 for v in rr],'o-',c=col,ms=3,label=lab)
zero(a);panel(a,'a','Old probe: a projection can vanish');a.set(xlabel=r'$U$ (meV)',ylabel=r'$C_{iX}$ [$10^{-5}$ (nm V T)$^{-1}$]');a.legend(fontsize=7)
rr=sorted([v for v in D['km_normal'] if v['t_meV']==10],key=lambda v:v['U_top_minus_bottom_meV'])
layers(b,rr,'U_top_minus_bottom_meV',1e-5);zero(b);panel(b,'b','New probe: wider detuning range');b.set(xlabel=r'$U$ (meV)',ylabel=unit5);b.legend(fontsize=7,loc='lower left')
fig.text(.5,.96,r'Fixed KM registry $\lambda=1/4$; $t_\perp=10$ meV',ha='center',fontsize=8)
save(fig,'FigS6_probe_and_range')
''' + script[end:]
script = script.replace("print('Eight figures exported (3 main, 5 supplementary).')", "print('Nine figures exported (3 main, 6 supplementary).')")
(R / 'plot_revision.py').write_text(script)

selections = {
    'fine_input_rows': len(D['km_normal']),
    'Fig1b': {'predicate': 'U=0', 'rows': sum(v['U_top_minus_bottom_meV']==0 for v in D['km_normal'])},
    'Fig1c': {'predicate': 't=10 and 0<=U<=5', 'rows': sum(v['t_meV']==10 and v['U_top_minus_bottom_meV']<=5 for v in D['km_normal']), 'reason': 'first total zero while both layer signs are retained; full scan in S6b'},
    'Fig1d': {'predicate': 'mass=45, lambda in {.24,.25,.26}, t=10, 3<=U<=5, common coarse grid', 'reason': 'zoom on nearby total zeros; all raw points archived'},
    'FigS6a': {'rows': 4, 'source': 'tensor_check.json; lambda=.25', 'reason': 'independent full in-plane tensor, native X probe; no full vector reconstructed from scalar data'},
    'statistics': 'deterministic calculations; no error bars or statistical sample claims',
}
(R / 'figure-selection.json').write_text(json.dumps(selections, indent=2) + '\n')
print('Figure source assembled; no response calculations run by this script.')
