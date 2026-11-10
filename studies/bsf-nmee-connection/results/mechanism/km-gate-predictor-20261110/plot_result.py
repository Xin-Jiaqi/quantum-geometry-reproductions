"""Physics figure: unequal cancellation gives layer-selective gate sensitivity.

Python quantitative grid, 183x82 mm, editable SVG/PDF and 600dpi PNG.
All retained source records in source_data.csv; no fits/statistical errors.
"""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'prediction-fine.json').read_text())['rows'][0]
rows=[x for x in json.loads((ROOT/'verification-fine.json').read_text())['rows'] if x['t_meV']==10]
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7})
fig,axs=plt.subplots(1,2,figsize=(183/25.4,82/25.4))
fig.subplots_adjust(left=.10,right=.98,bottom=.26,top=.83,wspace=.38)
colors=['#346995','#BC702C']
for ax,panel in zip(axs,'ab'):
 ax.axhline(0,lw=.7,color='#777777');ax.text(-.16,1.08,panel,transform=ax.transAxes,fontsize=10,fontweight='bold')
for i,(layer,c) in enumerate(zip(['Lower layer','Upper layer'],colors)):
 axs[0].bar(np.arange(3)+(i-.5)*.30,[1000*r[x][i] for x in ['energy_shift','matrix_change','tangent']],width=.28,color=c,label=layer)
 states=sorted(rows,key=lambda x:x['U_meV']);u=np.array([x['U_meV'] for x in states]);delta=np.array([x['chi'][i]-r['chi'][i] for x in states])
 axs[1].plot(u,u/1000*r['tangent'][i]*1e6,'--',color=c,lw=1)
 axs[1].plot(u,delta*1e6,'o',color=c,ms=4)
axs[0].set_xticks(np.arange(3),['Transition\nenergies','Spin/electric\nmatrix elements','Total'])
axs[0].set_ylabel(r'$\partial C_\ell/\partial U$'+'\n'+r'[$10^{-3}$ (nm V T eV)$^{-1}$]')
axs[0].set_title('Unequal cancellation within each layer',pad=14)
axs[0].legend(frameon=False,loc='lower right')
axs[1].set_xlabel('Upper-minus-lower energy offset U (meV)')
axs[1].set_ylabel(r'$C_\ell(U)-C_\ell(0)$'+'\n'+r'[$10^{-6}$ (nm V T)$^{-1}$]')
axs[1].set_title('Zero-bias prediction of gate response',pad=14)
axs[1].text(.5,-.36,'Dashed: predicted tangent    Points: finite-bias response',ha='center',transform=axs[1].transAxes,fontsize=7)
fig.text(.55,.96,r'Bilayer KM: $s_{\parallel,\ell}=C_\ell E_\parallel B_z$',ha='center',fontsize=10)
fig.savefig(ROOT/'km_gate_mechanism.png',dpi=600)
fig.savefig(ROOT/'km_gate_mechanism.svg')
fig.savefig(ROOT/'km_gate_mechanism.pdf')
plt.close(fig)
p=ROOT/'km_gate_mechanism.svg';p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
