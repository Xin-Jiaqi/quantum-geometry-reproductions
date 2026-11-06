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
 out=R/'figures-20261106';out.mkdir(exist_ok=True)
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
 u=np.array([-1/6,7/12])@lattice
 coords=np.array([np.array([i,j])@lattice+q for i in range(3) for j in range(3) for q in ([0,0],[0,1])])
 allxy=[]
 for layer,off,z in [('lower',u,0),('upper',np.zeros(2),3.3)]:
  p=coords+off;xy=np.column_stack((p[:,0]+.25*p[:,1],.31*p[:,1]+z));allxy.append(xy)
  for i in range(len(p)):
   for j in range(i):
    if abs(np.linalg.norm(p[i]-p[j])-1)<1e-6:ax.plot(xy[[i,j],0],xy[[i,j],1],c='.58',lw=.65,zorder=1)
  for parity,c in [(0,'#5D7AA8'),(1,'#E2AF58')]:ax.scatter(xy[parity::2,0],xy[parity::2,1],s=23,c=c,edgecolors='white',lw=.4,zorder=3)
 for i,j in [(4,4),(8,9),(12,12)]:
  ax.plot([allxy[0][i,0],allxy[1][j,0]],[allxy[0][i,1],allxy[1][j,1]],'--',c='#805891',lw=.8,zorder=0)
 ax.text(-2.5,5.7,r'Upper: $+\lambda_R,\ +U/2$',color=T,fontsize=8)
 ax.text(-2.5,-.5,r'Lower: $-\lambda_R,\ -U/2$',color=B,fontsize=8)
 ax.text(-2.2,2.9,r'$t_\perp T$',color='#805891',fontsize=9)
 ax.text(3.0,4.7,r'$t_0,\lambda_{\rm SO}$',fontsize=8)
 ax.text(-2.5,1.4,r'$a:\ +m$'+'\n'+r'$b:\ -m$',fontsize=8)
 ax.set(xlim=(-2.7,5),ylim=(-.9,6.1));ax.axis('off')

fig,axs=plt.subplots(2,2,figsize=(7,5.7));fig.subplots_adjust(left=.11,right=.97,bottom=.09,top=.91,wspace=.40,hspace=.60)
a,b,c,d=axs.flat;km_geometry(a);panel(a,'a','Mirror-oriented Kane–Mele bilayer')
r=sorted([v for v in D['km'] if v['U_top_minus_bottom_meV']==0],key=lambda v:v['t_meV'])
layers(b,r,'t_meV',1e-5);panel(b,'b','Interlayer hybridization');b.set(xlabel=r'$t_\perp$ (meV)',ylabel=unit5);b.legend(loc='upper left');b.text(.96,.05,r'$U=0$',transform=b.transAxes,ha='right')
r2=sorted([v for v in D['km'] if v['t_meV']==10],key=lambda v:v['U_top_minus_bottom_meV'])
layers(c,r2,'U_top_minus_bottom_meV',1e-5);panel(c,'c','Layer potential difference');c.set(xlabel=r'$U=\epsilon_{\rm upper}-\epsilon_{\rm lower}$ (meV)',ylabel=unit5);c.text(.03,.05,r'$t_\perp=10$ meV',transform=c.transAxes)
ref=next(v for v in r if v['t_meV']==10)
for key,col,lab in [('bottom_positive_k',T,r'$\Delta C_{+}$'),('bottom_negative_k',B,r'$\Delta C_{-}$'),('bottom',N,r'$\Delta C_{\rm lower}$')]:
 vals=[v for v in r if v.get(key)!='']
 d.plot([v['t_meV'] for v in vals],[(v[key]-ref[key])*1e5 for v in vals],'o-',c=col,ms=3,lw=1.1,label=lab)
zero(d);panel(d,'d','Redistribution within the lower layer');d.set(xlabel=r'$t_\perp$ (meV)',ylabel=r'$\Delta C_{xx}$ [same units as (b)]');d.legend(loc='upper left',ncol=1);d.text(.96,.06,r'Reference: $t_\perp=10$ meV',transform=d.transAxes,ha='right',fontsize=7)
save(fig,'Fig1_KM')

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

fig=plt.figure(figsize=(7,5.5));gs=fig.add_gridspec(2,6,left=.11,right=.97,bottom=.09,top=.91,hspace=.50,wspace=.72)
for i,dx in enumerate([0,1,2]):
 ax=fig.add_subplot(gs[0,2*i:2*i+2]);material(ax,dx);ax.set_title(r'$\Delta d='+str(dx)+r'\ \mathrm{\AA}$',fontsize=9,pad=0)
 if i==0:ax.text(-.15,1.1,'a',transform=ax.transAxes,fontweight='bold',fontsize=10);ax.legend(loc='lower center',ncol=3,fontsize=7,handletextpad=.2,columnspacing=.6)
ax=fig.add_subplot(gs[1,:3]);arr=np.array(D['distance']);ax.plot(arr[:,0],arr[:,1],'o-',c=B,ms=4,lw=1.2);zero(ax);panel(ax,'b','Lower-layer spin response');ax.set(xlabel=r'Added separation $\Delta d$ ($\mathrm{\AA}$)',ylabel=unit9)
ax=fig.add_subplot(gs[1,3:]);r=D['decomposition']
for i,v in enumerate(r):
 ax.bar(i-.17,v['deep'],width=.29,color=B,label='Deeper valence bands' if i==0 else None)
 ax.bar(i+.17,v['edge'],width=.29,color=T,label='Upper valence bands' if i==0 else None)
 ax.plot(i,v['total'],'D',c=N,ms=4,label='Lower-layer total' if i==0 else None)
zero(ax);panel(ax,'c','Occupied-band contributions');ax.set(xticks=[0,1],xticklabels=[r'$\Delta d=0$',r'$\Delta d=1\ \mathrm{\AA}$'],ylabel=unit9,ylim=(-36,35));ax.legend(loc='lower center',fontsize=6.5)
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
print('Six figures exported.')
