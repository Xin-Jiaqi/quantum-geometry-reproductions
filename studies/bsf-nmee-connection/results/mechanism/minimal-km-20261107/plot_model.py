from pathlib import Path
import json,csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
p=Path(__file__).parent
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False})
a=json.loads((p/'controlled.json').read_text())['rows'];b=json.loads((p/'main-thresholds.json').read_text())
series=[a['main_t'+str(g)] for g in [0.,.04,.08,.12,.16]]
rows=[]
for r in series:rows.append(['a','calculated',r['parameters']['textra']*1000,r['C'][0]*1e5,'controlled.json:'+r['parameters']['name']])
for r in b['training']:rows.append(['b','training',r['energy']*1000,r['t']*1000,'main-thresholds.json:training'])
for r in b['predictions']:rows.append(['b','held-out',r['energy']*1000,r['actual_t']*1000,'main-thresholds.json:predictions'])
with (p/'figure_source.csv').open('w',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['panel','series','x','y','source']);w.writerows(rows)
fig,ax=plt.subplots(1,2,figsize=(7.2047244,2.7559055),layout='constrained')
blue='#305D8A';orange='#C27A40'
ax[0].axhline(0,color='.6',lw=.7)
ax[0].plot([r['parameters']['textra']*1000 for r in series],[r['C'][0]*1e5 for r in series],'-o',color=blue,lw=1.2,ms=4)
ax[0].set(xlabel='Valence–empty interlayer coupling g (meV)',ylabel=r'$C_{\mathrm{lower}}$ [$10^{-5}$ (nm V T)$^{-1}$]')
ax[0].text(.04,.95,r'$\Delta=120$ meV, $\rho=0.1$',transform=ax[0].transAxes,va='top',fontsize=8)
ax[0].set_xlim(-5,170)
x=np.linspace(.12,.5,100)
ax[1].plot(x*1000,np.sqrt(b['scale']*x)*1000,'--',color='.45',lw=1,label=r'$g_c\propto\sqrt{\Delta}$ estimate')
ax[1].plot([r['energy']*1000 for r in b['training']],[r['t']*1000 for r in b['training']],'o',color=blue,ms=4,label='Reference calculations')
ax[1].plot([r['energy']*1000 for r in b['predictions']],[r['actual_t']*1000 for r in b['predictions']],'s',color=orange,ms=5,label='Predicted cases')
ax[1].set(xlabel=r'Additional empty-state energy $\Delta$ (meV)',ylabel=r'Reversal threshold $g_c$ (meV)')
ax[1].legend(loc='upper left',fontsize=7)
for label,a in zip(['a','b'],ax):a.text(-.14,1.03,label,transform=a.transAxes,weight='bold',fontsize=10)
fig.savefig(p/'minimal_model_mechanism.svg')
fig.savefig(p/'minimal_model_mechanism.pdf')
fig.savefig(p/'minimal_model_mechanism.png',dpi=600)

svg=p/'minimal_model_mechanism.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
