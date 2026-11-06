"""Model, material, and auxiliary figures from archived data; no new response simulation."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;D=json.loads((R/'figure_data.json').read_text())
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
B='#286C9C';T='#C87D36';N='#242B35';A='#75439E';P='#14867B'
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=12)
 ax.text(-.16,1.07,letter,transform=ax.transAxes,fontsize=10,fontweight='bold')
def zero(ax):ax.axhline(0,c='.65',lw=.7,zorder=0)
def save(fig,name):
 out=R/'figures-20261106b';out.mkdir(exist_ok=True)
 fig.savefig(out/f'{name}.png',dpi=320)
 fig.savefig(out/f'{name}.pdf')
 fig.savefig(out/f'{name}.svg')
 svg=out/f'{name}.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
 plt.close(fig)
def layers(ax,rows,xkey,scale):
 for key,c,l in [('bottom',B,'Lower'),('top',T,'Upper'),('total',N,'Total')]:
  ax.plot([v[xkey] for v in rows],[v[key]/scale for v in rows],'o-',c=c,ms=3,lw=1.2,label=l)
 zero(ax)
unit9=r'$C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]'
unit5=r'$C_{xx}$ [$10^{-5}$ (nm V T)$^{-1}$]'

def km_geometry(ax):
 lattice=np.array([[-np.sqrt(3)/2,1.5],[np.sqrt(3)/2,1.5]])
 base=np.array([np.array([i,j])@lattice+q for i in range(3) for j in range(3) for q in ([0,0],[0,1])])
 for panel_x,shift,label in [(0,np.array([-1/3,2/3])@lattice,'AB'),(7,np.array([1/3,1/3])@lattice,'BA')]:
  for name,off,marker,fc,ec in [('upper',np.zeros(2),'o',T,T),('lower',shift,'s','white',B)]:
   xy=base+off+np.array([panel_x,0])
   for i in range(len(xy)):
    for j in range(i):
     if abs(np.linalg.norm(xy[i]-xy[j])-1)<1e-6:ax.plot(xy[[i,j],0],xy[[i,j],1],c=ec,lw=.6,alpha=.45)
   ax.scatter(xy[:,0],xy[:,1],s=22,marker=marker,facecolors=fc,edgecolors=ec,lw=.8,zorder=3)
  ax.text(panel_x+1,8.3,label,ha='center',fontsize=9)
 ax.annotate('',xy=(5.4,2.6),xytext=(3.5,2.6),arrowprops={'arrowstyle':'<->','color':N,'lw':.8})
 ax.text(4.4,3.1,r'$m_{001}$',ha='center',fontsize=8)
 ax.text(-1.7,-1.2,r'Upper: $+\lambda_R,+U/2$',c=T,fontsize=7)
 ax.text(-1.7,-2.0,r'Lower: $-\lambda_R,-U/2$',c=B,fontsize=7)
 ax.set(xlim=(-2.5,11.5),ylim=(-2.4,9.0));ax.set_aspect('equal');ax.axis('off')

fig,axs=plt.subplots(2,2,figsize=(7,5.7));fig.subplots_adjust(left=.11,right=.97,bottom=.09,top=.91,wspace=.40,hspace=.60)
a,b,c,d=axs.flat;km_geometry(a);panel(a,'a','AB and its mirror partner')
r=sorted([v for v in D['km'] if v['U_top_minus_bottom_meV']==0],key=lambda v:v['t_meV'])
layers(b,r,'t_meV',1e-5);panel(b,'b','Hybridization at AB');b.set(xlabel=r'$t_\perp$ (meV)',ylabel=unit5);b.legend(loc='upper left',fontsize=7);b.text(.96,.05,r'$U=0$',transform=b.transAxes,ha='right')
r2=sorted([v for v in D['km'] if v['t_meV']==10],key=lambda v:v['U_top_minus_bottom_meV'])
layers(c,r2,'U_top_minus_bottom_meV',1e-5);panel(c,'c','Layer detuning at AB');c.set(xlabel=r'$U$ (meV)',ylabel=unit5);c.text(.03,.94,r'$t_\perp=10$ meV',transform=c.transAxes,va='top')
ar=next(v for v in D['km'] if v['t_meV']==10 and v['U_top_minus_bottom_meV']==0)
br=next(v for v in D['km_ba'] if v['t_meV']==10 and v['U_top_minus_bottom_meV']==0)
for j,(key,col,lab) in enumerate([('bottom',B,'Lower'),('top',T,'Upper'),('total',N,'Total')]):
 d.bar(np.array([0,1])+(j-1)*.22,[ar[key]*1e5,br[key]*1e5],width=.20,color=col,label=lab)
zero(d);panel(d,'d','Layer exchange under reflection');d.set(xticks=[0,1],xticklabels=['AB','BA'],ylabel=unit5);d.text(.97,.94,r'$t_\perp=10$ meV'+'\n'+r'$U=0$',transform=d.transAxes,ha='right',va='top',fontsize=7)
save(fig,'Fig1_KM')

# The previously displayed generic displacement remains available, unchanged, in SI.
fig,(b,c,d)=plt.subplots(1,3,figsize=(7,2.85));fig.subplots_adjust(left=.09,right=.98,bottom=.23,top=.82,wspace=.56)
r=sorted([v for v in D['km_quarter'] if v['U_top_minus_bottom_meV']==0],key=lambda v:v['t_meV'])
layers(b,r,'t_meV',1e-5);panel(b,'a','Hybridization');b.set(xlabel=r'$t_\perp$ (meV)',ylabel=unit5);b.legend(fontsize=6)
r2=sorted([v for v in D['km_quarter'] if v['t_meV']==10],key=lambda v:v['U_top_minus_bottom_meV'])
layers(c,r2,'U_top_minus_bottom_meV',1e-5);panel(c,'b','Layer detuning');c.set(xlabel=r'$U$ (meV)',ylabel=unit5)
ref=next(v for v in r if v['t_meV']==10)
for key,col,lab in [('bottom_positive_k',T,r'$\Delta C_{+}$'),('bottom_negative_k',B,r'$\Delta C_{-}$'),('bottom',N,r'$\Delta C_{\rm lower}$')]:
 vals=[v for v in r if v.get(key)!=''];d.plot([v['t_meV'] for v in vals],[(v[key]-ref[key])*1e5 for v in vals],'o-',c=col,ms=3,lw=1.1,label=lab)
zero(d);panel(d,'c','Signed k contributions');d.set(xlabel=r'$t_\perp$ (meV)',ylabel=r'$\Delta C_{xx}$ [same units]');d.legend(fontsize=6)
fig.text(.5,.99,'Fixed stacking: one quarter along AB–BA',ha='center',va='top',fontsize=8)
save(fig,'FigS4_generic_stacking')

def material(ax,extra=0,distort=False):
 cell=np.array(D['structure_A']['cell']);pos=np.array(D['structure_A']['positions']);typ=D['structure_A']['species'];pts=[];sp=[]
 mid=cell[2,2]/2
 for i in range(3):
  for j in range(2):
   for f,s in zip(pos,typ):
    p=(f+np.array([i,j,0]))@cell;p[2]-=mid
    p[2]+=extra/2 if p[2]>0 else -extra/2
    pts.append(p);sp.append(s)
 pts=np.array(pts);xy=np.column_stack((pts[:,0]+.18*pts[:,1],pts[:,2]+.13*pts[:,1]))
 for i in range(len(pts)):
  for j in range(i):
   if (sp[i]=='Mo')!=(sp[j]=='Mo') and np.linalg.norm(pts[i]-pts[j])<2.85:ax.plot(xy[[i,j],0],xy[[i,j],1],c='.64',lw=.7,zorder=0)
 for s,col,size in [('Mo','#7668A0',38),('S','#DDB545',27),('Se','#46A493',33)]:
  mask=np.array(sp)==s;ax.scatter(xy[mask,0],xy[mask,1],s=size,c=col,edgecolors='white',lw=.35,label=s,zorder=2)
 if distort:
  idx=next(i for i,s in enumerate(sp) if s=='Se' and pts[i,2]<0);x,y=xy[idx];ax.annotate('',xy=(x,y+.9),xytext=(x,y-.8),arrowprops={'arrowstyle':'->','color':'#BC4D47','lw':1.5});ax.text(x+1,y-1.5,'0.01 Å\n(exaggerated)',fontsize=7,color='#BC4D47')
 ax.set_aspect('equal');ax.axis('off');ax.set(xlim=(-2.3,10),ylim=(-6.5,7.4))

fig=plt.figure(figsize=(7,6));gs=fig.add_gridspec(2,2,left=.11,right=.97,bottom=.09,top=.92,hspace=.62,wspace=.40)
sub=gs[0,0].subgridspec(1,2,wspace=.05)
for i,dx in enumerate([0,1]):
 ax=fig.add_subplot(sub[0,i]);material(ax,dx);ax.set_title(r'$\Delta d='+str(dx)+r'\ \mathrm{\AA}$',fontsize=9,pad=0)
 if i==0:ax.text(-.22,1.12,'a',transform=ax.transAxes,fontweight='bold',fontsize=10)
 if i==1:ax.legend(loc='lower center',bbox_to_anchor=(.5,-.12),ncol=3,fontsize=6,handletextpad=.1,columnspacing=.3)
ax=fig.add_subplot(gs[0,1]);arr=np.array(D['distance']);ax.plot(arr[:,0],arr[:,1],'o-',c=B,ms=4,lw=1.2,label='Lower layer');ax.axhline(D['monolayer_reference'],c='.4',ls='--',lw=.9,label='Isolated monolayer');zero(ax);panel(ax,'b','Interface-induced reversal');ax.set(xlabel=r'Added separation $\Delta d$ ($\mathrm{\AA}$)',ylabel=unit9);ax.legend(fontsize=6,loc='lower right')
ax=fig.add_subplot(gs[1,0]);r=D['decomposition']
for i,v in enumerate(r):
 ax.bar(i-.16,v['deep'],width=.29,color=B,label=r'$n=1$–$36$' if i==0 else None)
 ax.bar(i+.16,v['edge'],width=.29,color=T,label=r'$n=37$–$48$' if i==0 else None)
 ax.plot(i,v['total'],'D',c=N,ms=4,label='Lower-layer total' if i==0 else None)
zero(ax);panel(ax,'c','Occupied-band contributions');ax.set(xticks=[0,1],xticklabels=['0','1'],xlabel=r'$\Delta d$ ($\mathrm{\AA}$)',ylabel=unit9);ax.legend(fontsize=6,loc='upper left')
ax=fig.add_subplot(gs[1,1])
for i,case in enumerate(D['bands']):
 for j,p in enumerate(case['points']):
  x=3*i+j
  for n,e in enumerate(p['energy_relative_vbm'],1):
   if n<21:continue
   col=B if n<=36 else T if n<=48 else '.65'
   ax.hlines(e,x-.22,x+.22,color=col,lw=1.1)
ax.axhline(0,c='.85',lw=.7);ax.axvline(2,c='.85',lw=.7)
panel(ax,'d','Energy levels of the band groups');ax.set(xticks=[0,1,3,4],xticklabels=[r'$\Gamma$',r'$K$',r'$\Gamma$',r'$K$'],ylabel=r'$\varepsilon-E_{\rm VBM}$ (eV)',ylim=(-8,4),xlim=(-.6,4.6))
ax.text(.21,-.23,r'$\Delta d=0$',transform=ax.transAxes,ha='center',fontsize=7);ax.text(.79,-.23,r'$\Delta d=1$ Å',transform=ax.transAxes,ha='center',fontsize=7)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],c=B,label='21–36'),Line2D([0],[0],c=T,label='37–48'),Line2D([0],[0],c='.65',label='49–56')],fontsize=6,ncol=3,loc='upper center',handlelength=.9,columnspacing=.7)
save(fig,'Fig2_MoSSe_distance')

fig,axs=plt.subplots(1,3,figsize=(7,2.85));fig.subplots_adjust(left=.085,right=.985,bottom=.20,top=.79,wspace=.65)
a,b,c=axs
for key,lab,color in [('projected44','Mo d–S/Se p',B),('embedded112','+ remaining bands',T)]:
 rows=sorted(D[key],key=lambda v:v.get('hopping',v.get('t')));a.plot([v.get('hopping',v.get('t')) for v in rows],[v['C'] for v in rows],'o-',c=color,ms=3,lw=1,label=lab)
zero(a);panel(a,'a','Hybridization only');a.set(xlabel=r'Interlayer scale $\alpha$',ylabel=r'Lower-layer '+unit9);a.legend(fontsize=6,loc='upper right')
pa=D['mirror']['paired'];fr=[v['F_A'] for v in pa]
for idx,col,lab in [(1,B,'Lower'),(2,T,'Upper'),(0,N,'Total')]:b.plot(fr,[v['A'][idx] for v in pa],'o-',c=col,ms=3,lw=1,label=lab)
zero(b);panel(b,'b','Vertical electric field');b.set(xlabel=r'$F_z$ (V/nm)',ylabel=unit9);b.legend(loc='upper left',bbox_to_anchor=(.02,.73),fontsize=6)
for dom,col,sgn in [('A',A,1),('B',P,-1)]:
 pairs=sorted([(sgn*v['F_A'],v[dom][0]) for v in pa]);c.plot(*np.array(pairs).T,'o-',c=col,ms=3,lw=1,label='Stacking '+dom)
zero(c);panel(c,'c','Mirror stacking partners');c.set(xlabel=r'$F_z$ (V/nm)',ylabel=r'Total '+unit9);c.legend(loc='lower right',fontsize=6)
save(fig,'Fig3_MoSSe_controls')

fig,(a,b)=plt.subplots(1,2,figsize=(7,3.05));fig.subplots_adjust(left=.08,right=.97,bottom=.18,top=.80,wspace=.35)
material(a,0,True);panel(a,'a','Single-sided structural perturbation')
for dom,col in [('A',A),('B',P)]:
 vals=D['device']['lower_Se_inward'][dom]['rows'];b.plot([r[0] for r in vals],[r[1][0] for r in vals],'o-',c=col,ms=4,lw=1.1,label=dom+'; displaced Se')
 vals=sorted([(v['F_A']*(1 if dom=='A' else -1),v[dom][0]) for v in pa]);b.plot(*np.array(vals).T,'--',c=col,alpha=.65,lw=1,label=dom+'; pristine')
zero(b);panel(b,'b','Shift of the response zero');b.set(xlabel=r'$F_z$ (V/nm)',ylabel=r'Total '+unit9);b.legend(fontsize=6,loc='lower right')
save(fig,'FigS1_environment')

fig,(a,b)=plt.subplots(1,2,figsize=(7,2.95));fig.subplots_adjust(left=.11,right=.97,bottom=.20,top=.80,wspace=.40)
for ax,mode,title,letter in [(a,'above_aperture','Probe above the edge','a'),(b,'aperture','Laterally displaced probe','b')]:
 for dom,col in [('A',A),('B',P)]:
  vals=sorted([v for v in D['edges'] if v['mode']==mode and v['domain']==dom],key=lambda v:v['F']);ax.plot([v['F'] for v in vals],[v['field_aT'] for v in vals],'o-',c=col,ms=3,lw=1.1,label='Stacking '+dom)
 zero(ax);panel(ax,letter,title);ax.set(xlabel=r'$F_z$ (V/nm)',ylabel=r'Calculated stray field (aT)');ax.legend(fontsize=7)
fig.text(.5,.97,r'Ideal spin-only edge model: $E_x=1$ V/μm, $B_z=1$ T; probe height / width = 50 nm',ha='center',va='top',fontsize=8)
save(fig,'FigS2_probe_weight')

fig,(a,b)=plt.subplots(1,2,figsize=(7,2.95));fig.subplots_adjust(left=.11,right=.97,bottom=.20,top=.80,wspace=.42)
for dom,col in [('A',A),('B',P)]:
 vals=sorted([(v['F_A']*(1 if dom=='A' else -1),v[dom][0]) for v in pa]);a.plot(*np.array(vals).T,'o-',c=col,ms=3,label=dom+'; original probe')
 vals=D['device']['weak'][dom]['rows'];a.plot([v[0] for v in vals],[v[1][0] for v in vals],'s--',mfc='white',mec=col,c=col,ms=4,label=dom+'; half probe')
zero(a);panel(a,'a','Smaller probing fields');a.set(xlabel=r'$F_z$ (V/nm)',ylabel=r'Total '+unit9);a.legend(fontsize=6,loc='upper left')
x=[0,1];y=[next(v[1] for v in D['distance'] if v[0]==q) for q in x];b.plot(x,y,'o-',c=B,ms=5,label=r'$6\times6$');b.plot(x,[v['total'] for v in D['decomposition']],'s--',c=T,mfc='white',ms=5,label=r'$9\times9$')
zero(b);panel(b,'b','Endpoint signs on two k meshes');b.set(xlabel=r'$\Delta d$ ($\mathrm{\AA}$)',ylabel=r'Lower-layer '+unit9,xticks=[0,1]);b.legend(loc='upper left')
save(fig,'FigS3_numerical_checks')


# Orbital character only at archived high-symmetry points; no invented band paths.
fig,axs=plt.subplots(2,3,figsize=(7,5.1));fig.subplots_adjust(left=.09,right=.98,bottom=.13,top=.87,wspace=.22,hspace=.36)
channels=[('Mo d','#7668A0',lambda k:'-Mo-d' in k),(r'S $p_z$','#C49D2D',lambda k:'-S-pz' in k),(r'Se $p$','#238A70',lambda k:'-Se-p' in k)]
for row,case in enumerate(D['orbital_levels']):
 pts={p['kname']:p for p in case['points']}
 for col,(label,color,test) in enumerate(channels):
  ax=axs[row,col]
  for x,k in enumerate(['G','M','K','Kp']):
   p=pts[k];e=np.array(p['energy'])-case['vbm'];den=np.sum(list(p['projected'].values()),axis=0);w=np.sum([v for key,v in p['projected'].items() if test(key)],axis=0)/den
   ax.scatter(np.full(len(e),x),e,s=3,c='.8',zorder=0)
   ax.scatter(np.full(len(e),x),e,s=110*w,facecolors='none',edgecolors=color,lw=.7)
  ax.axhline(0,c='.8',lw=.6);ax.set(xticks=range(4),xticklabels=[r'$\Gamma$',r'$M$',r'$K$',r'$K^\prime$'],ylim=(-5.3,4.5),xlim=(-.6,3.6))
  if row==0:ax.set_title(label,fontsize=9)
  if col==0:ax.set_ylabel(r'$\varepsilon-E_{\rm VBM}^{*}$ (eV)');ax.text(-.30,.5,r'$\Delta d='+str(case['distance_A'])+r'$ Å',transform=ax.transAxes,rotation=90,ha='center',va='center',fontsize=8)
  else:ax.tick_params(labelleft=False)
  ax.text(-.12,1.03,chr(97+row*3+col),transform=ax.transAxes,fontweight='bold',fontsize=9)
fig.text(.5,.99,'Orbital character at sampled k points (bands 39–56)',ha='center',va='top',fontsize=9)
fig.text(.5,.02,'Circle area: fraction of the summed PAW-sphere weight; both layers included',ha='center',fontsize=7)
save(fig,'FigS5_orbital_levels')
print('Eight figures exported (3 main, 5 supplementary).')
