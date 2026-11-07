"""Plot archived static-spin responses; no model fitting or simulated observations."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
ROWS=list(csv.DictReader((ROOT/'source_data.csv').open()))
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'legend.frameon':False,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
COL={'bottom':'#286C9C','top':'#C87D36','total':'#252B35'}
NAMES={'bottom':'Lower layer','top':'Upper layer','total':'Net'}
def data(group,layer):
 r=sorted((v for v in ROWS if v['group']==group and v['layer']==layer),key=lambda v:float(v['x']))
 return np.array([float(v['x']) for v in r]),np.array([float(v['C']) for v in r])
def curves(ax,group,scale,layers=('bottom','top','total')):
 for layer in layers:
  x,y=data(group,layer)
  ax.plot(x,y/scale,'o-',c=COL[layer],ms=3,lw=1.25,label=NAMES[layer])
 ax.axhline(0,c='.6',lw=.7,zorder=0)
def panel(ax,label,title):
 ax.text(-.16,1.06,label,transform=ax.transAxes,fontweight='bold',fontsize=10)
 ax.set_title(title,loc='left',pad=10)
def export(fig,name):
 fig.savefig(ROOT/f'{name}.png',dpi=320)
 fig.savefig(ROOT/f'{name}.pdf')
 fig.savefig(ROOT/f'{name}.svg')
 fig.savefig(ROOT/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(7,5.45))
fig.subplots_adjust(left=.11,right=.97,bottom=.09,top=.86,wspace=.38,hspace=.56)
a,b,c,d=axs.flat
curves(a,'KM_t',1e-5)
a.set(xlabel=r'Interlayer hybridization $t_\perp$ (meV)',ylabel=r'$C_{xx}$ [$10^{-5}$ (nm V T)$^{-1}$]')
panel(a,'a','KM: a layer reverses')
a.legend(loc='upper left',ncol=1)
a.plot([10],[-2.9601943898546865],marker='o',ms=6,mfc='none',mec=COL['bottom'])
a.text(.98,.07,r'$U=0$',transform=a.transAxes,ha='right',fontsize=8)
curves(b,'KM_U',1e-5)
b.set(xlabel=r'Layer energy offset $U$ (meV)',ylabel=r'$C_{xx}$ [$10^{-5}$ (nm V T)$^{-1}$]')
panel(b,'b','KM: only the net response reverses')
b.text(.03,.06,r'$t_\perp=10$ meV',transform=b.transAxes,fontsize=8)
curves(c,'MoSSe_d',1e-9,layers=('bottom',))
c.set(xlabel=r'Added layer separation $\Delta d$ ($\mathrm{\AA}$)',ylabel=r'$C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]')
panel(c,'c','MoSSe: the same lower layer reverses')
c.text(.97,.08,'Fixed registry; zero control field',transform=c.transAxes,ha='right',fontsize=7)
curves(d,'MoSSe_A',1e-9)
d.set(xlabel=r'Vertical control field $F_z$ (V/nm)',ylabel=r'$C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]')
panel(d,'d','MoSSe: layer directions are retained')
d.axvline(.2441,color='.65',lw=.7,ls='--',zorder=0)
d.annotate(r'$F_0\simeq+0.24$',xy=(.2441,0),xytext=(.29,-2.5),fontsize=7,
 arrowprops={'arrowstyle':'-','color':'.4','lw':.7})
fig.text(.5,.98,r'Fixed probe directions: $\delta s_x^{\ell}=C_{xx}^{\ell}E_x B_z$',ha='center',va='top',fontsize=10)
export(fig,'layer_reversal_and_cancellation')

fig,(a,b)=plt.subplots(1,2,figsize=(7,3.25),gridspec_kw={'width_ratios':[1.15,1]})
fig.subplots_adjust(left=.04,right=.98,bottom=.19,top=.82,wspace=.36)
a.set(xlim=(0,10),ylim=(0,10));a.axis('off')
a.text(0,10.8,'a',fontweight='bold',fontsize=10)
a.text(5,11.4,r'Layer exchange + in-plane spin reversal',ha='center',fontsize=9)
for center,dom,shift in [(2.2,'A',r'$(1/3,\,2/3)$'),(7.7,'B',r'$(2/3,\,1/3)$')]:
 a.text(center,9.3,f'Stacking {dom}',ha='center',fontweight='bold',fontsize=9)
 a.text(center,8.3,shift,ha='center',fontsize=8)
 for layer,y,chem in [('top',6.6,'S–Mo–Se'),('bottom',2.6,'Se–Mo–S')]:
  x,vals=data('MoSSe_'+dom,layer);v=vals[np.argmin(abs(x))]/1e-9
  a.plot([center-1.75,center+1.75],[y-.38,y-.38],color=COL[layer],alpha=.22,lw=8,solid_capstyle='butt')
  length=abs(v)*.44
  start=center-length/2 if v>0 else center+length/2
  end=center+length/2 if v>0 else center-length/2
  a.annotate('',xy=(end,y),xytext=(start,y),arrowprops={'arrowstyle':'-|>','color':COL[layer],'lw':2,'mutation_scale':11})
  a.text(center,y+.6,f'{v:+.2f}',ha='center',color=COL[layer],fontsize=8)
  a.text(center,y-1.15,chem,ha='center',fontsize=7)
 a.text(center,.05,'Net '+('−0.416' if dom=='A' else '+0.416'),ha='center',fontsize=8)
a.text(4.95,4.8,r'$m_{001}$',ha='center',fontsize=10)
a.text(4.95,3.9,r'$z\to-z$',ha='center',fontsize=8)
a.text(4.95,-1.45,r'$E_x\to E_x,\quad B_z\to B_z,\quad s_x\to-s_x$',ha='center',fontsize=9)
for dom,color,marker in [('A','#7444A1','o'),('B','#198477','s')]:
 x,y=data('MoSSe_'+dom,'total');b.plot(x,y/1e-9,marker+'-',color=color,ms=4,lw=1.4,label=f'Stacking {dom}')
b.axhline(0,color='.5',lw=.7)
for x,lab in [(-.2440,r'$-0.24$'),(.2441,r'$+0.24$')]:
 b.axvline(x,color='.7',lw=.7,ls='--');b.text(x,-.38,lab,ha='center',fontsize=8)
b.set(xlabel=r'Vertical control field $F_z$ (V/nm)',ylabel=r'Net $C_{xx}$ [$10^{-9}$ (nm V T)$^{-1}$]',ylim=(-.49,.58))
b.text(-.20,1.15,'b',fontweight='bold',fontsize=10,transform=b.transAxes)
b.set_title('Opposite zero-crossing fields',loc='left',pad=16)
b.legend(loc='upper left',fontsize=7)
fig.text(.035,.025,r'Layer values: $10^{-9}$ (nm V T)$^{-1}$; arrows show induced spin for positive $E_xB_z$.',fontsize=7)
export(fig,'mirror_domains')
print('Exported two figures in PNG/PDF/SVG/TIFF')
