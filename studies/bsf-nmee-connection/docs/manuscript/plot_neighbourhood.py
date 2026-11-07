"""Compare both NMEE reversal routes using all existing neighbourhood points."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'
OUT = ROOT/'附件'
OUT.mkdir(exist_ok=True)
with (DATA/'km_neighbourhood.csv').open() as stream:
    rows = [{k: (v if k in ['grid', 'source'] else float(v)) for k, v in row.items()}
            for row in csv.DictReader(stream)]
assert len(rows) == 37
assert len({r['grid'] for r in rows}) == 1
assert all(r['common_gap_meV'] > 0 for r in rows)
assert max(abs(r['bottom']+r['top']-r['total']) for r in rows) < 1e-12

plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,
    'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
    'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none',
})
colours = ['#75439E', '#242B35', '#14867B']
markers = ['o','s','^']
fig, axes = plt.subplots(2,3,figsize=(7.2,5.9))
fig.subplots_adjust(left=.095,right=.98,bottom=.10,top=.88,wspace=.47,hspace=.82)
fig.text(.095,.975,r'Nearby KM registries: $m=45$ meV',fontsize=9,va='top')
fig.text(.095,.465,r'Nearby onsite energies: $\lambda=0.25$',fontsize=9,va='top')

record = {}
for row_index, (key, values) in enumerate([
    ('lambda',[.24,.25,.26]), ('mass_meV',[44.9,45,45.1])
]):
    a,b,c = axes[row_index]
    for i, value in enumerate(values):
        if key == 'lambda':
            selected = [r for r in rows if r['lambda']==value and r['mass_meV']==45]
            label = rf'$\lambda={value:.2f}$'
        else:
            selected = [r for r in rows if r['mass_meV']==value and r['lambda']==.25]
            label = rf'$m={value:g}$ meV'
        hybrid = sorted([r for r in selected if r['U_top_minus_bottom_meV']==0],key=lambda r:r['t_meV'])
        gate = sorted([r for r in selected if r['t_meV']==10],key=lambda r:r['U_top_minus_bottom_meV'])
        assert hybrid[0]['bottom'] < 0 < hybrid[-1]['bottom']
        assert all(r['bottom'] < 0 < r['top'] for r in gate)
        assert gate[0]['total'] > 0 > gate[-1]['total']
        colour,marker = colours[i],markers[i]
        a.plot([r['t_meV'] for r in hybrid],[r['bottom']*1e5 for r in hybrid],
               marker+'-',color=colour,ms=3.5,lw=1.2,label=label)
        x = [r['U_top_minus_bottom_meV'] for r in gate]
        for component,style,fill in [('bottom','-','full'),('top','--','none')]:
            b.plot(x,[r[component]*1e5 for r in gate],linestyle=style,marker=marker,
                   fillstyle=fill,color=colour,ms=3,lw=1.1)
        near_zero = [r for r in gate if r['U_top_minus_bottom_meV'] >= 4]
        c.plot([r['U_top_minus_bottom_meV'] for r in near_zero],
               [r['total']*1e5 for r in near_zero],marker+'-',color=colour,ms=3.5,lw=1.2)
        record[f'{key}={value}'] = {'hybridization_points':len(hybrid),'gate_points':len(gate),
                                   'individual_gate_layer_signs':'lower negative; upper positive'}
    for column,ax in enumerate([a,b,c]):
        ax.axhline(0,color='.65',lw=.75,zorder=0)
        ax.text(-.19,1.10,chr(97+row_index*3+column),transform=ax.transAxes,
                fontsize=10,fontweight='bold')
        ax.set_title(['One layer reverses','Layers stay opposite','Total spin crosses zero'][column],loc='left',pad=12)
        ax.set_xlabel(r'$t_\perp$ (meV)' if column==0 else r'$U$ (meV)')
    a.set_ylabel(r'Lower-layer $C_{xx}$')
    b.set_ylabel(r'Layer $C_{xx}$')
    c.set_ylabel(r'Total $C_{xx}$')
    a.set(xticks=[10,15,20,25],ylim=(-4.6,4.4))
    b.set(xlim=(-.15,5.15),ylim=(-4.6,14))
    c.set(xlim=(3.95,5.05),ylim=(-1.2,2.1),xticks=[4,4.5,5])
    a.legend(loc='lower right',handlelength=1.6,labelspacing=.35)
    b.legend(handles=[Line2D([],[],color='.35',ls='-',label='Lower'),
                      Line2D([],[],color='.35',ls='--',label='Upper')],
             loc='center left',handlelength=1.6,labelspacing=.35)
    a.text(.03,.95,r'$U=0$',transform=a.transAxes,va='top',fontsize=7)
    c.text(.97,.94,r'$t_\perp=10$ meV',transform=c.transAxes,ha='right',va='top',fontsize=7)

fig.text(.5,.02,r'All $C_{xx}$ in $10^{-5}$ (nm V T)$^{-1}$; symbols are calculated points',ha='center',fontsize=8)
fig.savefig(OUT/'Fig02_KM_neighbourhood.png',dpi=320)
fig.savefig(OUT/'Fig02_KM_neighbourhood.pdf')
fig.savefig(OUT/'Fig02_KM_neighbourhood.svg')
svg=OUT/'Fig02_KM_neighbourhood.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
plt.close(fig)
(ROOT/'figure_mapping.json').write_text(json.dumps({'source_rows':len(rows), 'panels':record,
     'filter':'all supplied points used in applicable sweeps; total panels zoom to U>=4, all gate points are shown in layer panels',
     'shared_reference':'lambda=.25,m=45 shown in both rows',
     'interpolation':'straight display segments only; no interpolated simulation data'},indent=2)+'\n')
print('Combined six-panel comparison exported; every supplied point retained.')

# The existing directional control, labelled by physical probe directions.
tensors = [r for r in json.loads((DATA/'tensor_check.json').read_text()) if r['lambda']==.25]
tensors.sort(key=lambda r:r['U_meV'])
gate_rows=[]
for name in ['normal_main.csv','gate_extra_main.csv']:
    with (DATA/name).open() as stream:
        gate_rows.extend(r for r in csv.DictReader(stream) if float(r['t_meV'])==10)
gate_rows.sort(key=lambda r:float(r['U_top_minus_bottom_meV']))
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,3.05))
fig.subplots_adjust(left=.10,right=.97,bottom=.20,top=.80,wspace=.40)
for i,colour,label in [(0,'#242B35',r'$C_{XX}^{\rm tot}$'),(1,'#75439E',r'$C_{YX}^{\rm tot}$')]:
    a.plot([r['U_meV'] for r in tensors],[r['tensor']['total'][i][0]*1e5 for r in tensors],
           'o-',color=colour,ms=3,label=label)
for key,colour,label in [('bottom','#286C9C','Lower'),('top','#C87D36','Upper'),('total','#242B35','Total')]:
    b.plot([float(r['U_top_minus_bottom_meV']) for r in gate_rows],
           [float(r[key])*1e5 for r in gate_rows],'o-',color=colour,ms=3,label=label)
for letter,ax,title,ylabel in [
    ('a',a,r'$E_X$ along native $X$',r'$C_{iX}$ [$10^{-5}$ (nm V T)$^{-1}$]'),
    ('b',b,'Electric field along mirror normal',r'$C_{xx}$ [$10^{-5}$ (nm V T)$^{-1}$]')
]:
    ax.axhline(0,color='.65',lw=.7)
    ax.set(xlabel=r'$U$ (meV)',ylabel=ylabel)
    ax.set_title(title,loc='left',pad=12)
    ax.text(-.16,1.09,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
    ax.legend(fontsize=7,loc='best')
fig.text(.5,.96,r'Fixed KM registry $\lambda=1/4$; $t_\perp=10$ meV',ha='center',fontsize=8)
fig.savefig(OUT/'Fig04_KM_probe_direction.png',dpi=320)
fig.savefig(OUT/'Fig04_KM_probe_direction.pdf')
fig.savefig(OUT/'Fig04_KM_probe_direction.svg')
svg=OUT/'Fig04_KM_probe_direction.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
plt.close(fig)
print('Directional control relabelled by field orientation; response data unchanged.')
