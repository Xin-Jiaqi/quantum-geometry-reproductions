"""Exploratory, unfitted twelve-state KM extension and causal controls."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np

STUDY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(STUDY), str(STUDY.parents[1]/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditBilayer, integration_grid
from bsf_nmee.remote_empty_km import ingredients, selective_extra, assemble, flatten, response


def cases():
    items = []
    def add(name, **kwargs):
        item = dict(name=name, energy=.12, tlow=.01, textra=.12, rho=.1,
                    mass=.045, fraction=.25, site=0, rashba=.020,
                    mode='selective', flat=False)
        item.update(kwargs); items.append(item)
    for t in [0., .04, .08, .12, .16]:
        add(f'main_t{t}', textra=t)
    for rho in [0., .25, .5, 1.]:
        add(f'rho{rho}', rho=rho)
    for energy in [.25, .5]:
        add(f'energy{energy}', energy=energy)
    add('site_B', site=1)
    add('no_Rashba', rashba=0.)
    for channel in ['VC', 'CC', 'VV']:
        add('all_paths_'+channel, mode=channel)
    for t in [0., .12]:
        add('flat_'+str(t), textra=t, flat=True)
    for name, variation in [('mass44', dict(mass=.044)), ('mass46', dict(mass=.046)),
                            ('fraction20', dict(fraction=.20)), ('fraction30', dict(fraction=.30))]:
        for t in [0., .12]:
            add(name+'_'+str(t), textra=t, **variation)
    return items


def run(grid, step, selected=None):
    k, weights, _ = integration_grid(AuditBilayer(), *grid)
    items = [c for c in cases() if selected is None or c['name'] in selected]
    rows = {c['name']: dict(parameters=c, C=np.zeros(3), vmax=-np.inf, cmin=np.inf) for c in items}
    groups = sorted({(c['mass'],c['fraction'],c['site'],c['rashba']) for c in items})
    for mass, fraction, site, rashba in groups:
        group = [c for c in items if (c['mass'],c['fraction'],c['site'],c['rashba']) == (mass,fraction,site,rashba)]
        for start in range(0, len(k), 512):
            data = ingredients(k[start:start+512], mass=mass, fraction=fraction, site=site, rashba=rashba)
            for c in group:
                if c['mode']=='selective':
                    h,v=selective_extra(data,c['energy'],c['tlow'],c['textra'],c['rho'])
                else:
                    h,v=assemble(data,c['energy'],c['tlow'],c['textra'],c['mode'])
                if c['flat']:h,v=flatten(h,v)
                val,ev,ec=response(h,v,data[-3],data[-2],weights[start:start+512],step,data[-1])
                r=rows[c['name']];r['C']+=val
                r['vmax'],r['cmin']=max(r['vmax'],ev),min(r['cmin'],ec)
        print('Completed',mass,fraction,site,rashba,flush=True)
    for row in rows.values():
        row['gap_eV']=row.pop('cmin')-row.pop('vmax')
        row['insulating']=bool(row['gap_eV']>0)
        row['layer_sum_error']=float(row['C'][0]+row['C'][1]-row['C'][2])
        row['C']=row['C'].tolist()
    return dict(grid=grid,nk=len(k),step_T=step,units='(nm V T)^-1',rows=rows,
                scope='Unfitted KM plus one spinful empty orbital per layer; projected orbital-selective coupling; not a MoSSe downfolding or strain prediction.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--grid',nargs=3,type=int,default=[41,40,64])
    ap.add_argument('--step',type=float,default=.002);ap.add_argument('--select',nargs='+')
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    d=run(a.grid,a.step,a.select);a.out.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({n:[r['C'][0]*1e5,r['gap_eV']] for n,r in d['rows'].items()},indent=2))
