"""Plot only stored computed data; see docs/STAGE8_FIGURE_CONTRACT.md."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
    'legend.fontsize':7,'svg.fonttype':'none','pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.2})
COLORS=['#0072B2','#D55E00','#009E73','#6F4C9B']


def main():
    data=json.loads((ROOT/'results/bilayer_stage8_audit.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(7.204724,5.236220),layout='constrained')
    exported=[]
    def line(ax,panel,name,x,y,**kwargs):
        ax.plot(x,y,label=name,**kwargs)
        exported.extend({'panel':panel,'series':name,'x':float(a),'y':float(b)}
                        for a,b in zip(x,y))
    ax=axes[0,0]
    small=[r for r in data['sliding'] if 0 < r['u'] <= .2]
    ref=next(r for r in small if r['u']==.1)
    for key,label,color,marker in zip(
        ['dipole_over_e_d_per_cell','trace','difference'],
        [r'$P_z$',r'$C_{xx}+C_{yy}$',r'$C_{xx}-C_{yy}$'],COLORS,['o','s','^']):
        values=np.abs([r[key]/ref[key] for r in small])
        assert np.all(values>0)
        line(ax,'a',label,[r['u'] for r in small],values,color=color,
             marker=marker,fillstyle='none',markersize=4)
    ax.set(xscale='log',yscale='log',xlabel=r'Slide amplitude $u$',
           ylabel=r'$|Q(u)/Q(0.1)|$',title='Local powers, separately normalized')
    ax.legend(frameon=False)
    ax=axes[0,1]
    rows=[r for r in data['sliding'] if r['u']<=.4];x=[r['u'] for r in rows]
    series=[('Cxx',[r['C'][0][0] for r in rows],COLORS[0],'o','-'),
            ('Cyy',[r['C'][1][1] for r in rows],COLORS[1],'s','-'),
            ('trace / 2',[r['trace']/2 for r in rows],COLORS[2],None,'--'),
            ('difference / 2',[r['difference']/2 for r in rows],COLORS[3],None,':')]
    for name,y,color,marker,style in series:
        line(ax,'b',name,x,np.array(y)*1e6,color=color,marker=marker,
             markersize=3,linestyle=style)
    ax.axhline(0,color='.7',lw=.6)
    ax.set(xlabel=r'Slide amplitude $u$',ylabel=r'$C$ [$10^{-6}$ / (nm V T)]',
           title='Net spin response: two competing parts')
    ax.legend(frameon=False)
    ax=axes[1,0];rows=data['channels'];x=np.array([r['t_eV']*1e3 for r in rows])
    for key,name,color,marker in zip(['full','same_only','cross_only','relative_sign'],
        ['Full','Same only','Cross only','Relative sign reversed'],COLORS,['o','s','^','D']):
        y=np.array([r[key]['difference'] for r in rows])*1e7/x**2
        line(ax,'c',name,x,y,color=color,marker=marker,markersize=3,fillstyle='none')
    ax.axhline(0,color='.7',lw=.6)
    ax.set(xlabel=r'Interlayer hopping $t_\perp$ [meV]',
        ylabel=r'$(C_{xx}-C_{yy})/t_\perp^2$ [$10^{-7}$ / (nm V T meV$^2$)]',
        title='Two hopping paths control anisotropy')
    ax.legend(frameon=False)
    ax=axes[1,1];rows=[r for r in data['flat_spectrum'] if r['t_eV']==.005]
    for key,name,color,marker in [('physical','Physical spectrum',COLORS[0],'o'),
                                 ('flat','Fixed spectrum',COLORS[1],'s')]:
        line(ax,'d',name,[r['u'] for r in rows],
             [r[key]['C'][0][0]*1e6 for r in rows],color=color,marker=marker,markersize=4)
    ax.axhline(0,color='.7',lw=.6)
    ax.set(xlabel=r'Slide amplitude $u$',ylabel=r'$C_{xx}$ [$10^{-6}$ / (nm V T)]',
           title='Changing states matters; energy weights too')
    ax.legend(frameon=False)
    for label,ax in zip('abcd',axes.flat):
        ax.text(-.15,1.07,label,transform=ax.transAxes,fontweight='bold',fontsize=9)
    folder=ROOT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/'bilayer_stage8_physics.png',dpi=600)
    fig.savefig(folder/'bilayer_stage8_physics.svg')
    fig.savefig(folder/'bilayer_stage8_physics.pdf')
    plt.close(fig)
    with (ROOT/'results/bilayer_stage8_figure_data.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['panel','series','x','y'])
        writer.writeheader();writer.writerows(exported)


if __name__=='__main__':
    main()
