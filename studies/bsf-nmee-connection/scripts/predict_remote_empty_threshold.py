"""Freeze a leading virtual-mixing scale, then test unseen empty-state energies.

This is a conditional model scaling check, not a universal sign theorem or
an independent material prediction. The three training roots are not fitted
to MoSSe. The tested energies are not used to set the scale.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from test_remote_empty_mechanism import AuditBilayer, integration_grid, ingredients, selective_extra, response


def run(output, grid):
    k,w,_=integration_grid(AuditBilayer(),*grid)
    data=ingredients(k);rows=[]
    def calc(energy,t):
        h,v=selective_extra(data,energy,.01,t,.1)
        c,ev,ec=response(h,v,data[-3],data[-2],w,g=data[-1])
        if ec<=ev:raise ValueError('Insulating prediction crossed into a metal')
        rows.append(dict(energy=energy,t=t,C=c.tolist(),gap=ec-ev))
        return c[0]
    roots=[]
    for energy in [.12,.25,.5]:
        t=brentq(lambda x:calc(energy,x),.06,.25,xtol=1e-5)
        roots.append(dict(energy=energy,t=t,g2_over_Delta=t*t/energy))
        print('Training',roots[-1],flush=True)
    scale=float(np.median([r['g2_over_Delta'] for r in roots]))
    predictions=[dict(energy=e,predicted_t=np.sqrt(scale*e)) for e in [.18,.35]]
    record=dict(grid=grid,training=roots,scale=scale,predictions=predictions,
                parameters=dict(tlow=.01,rho=.1,mass=.045,fraction=.25,site=0))
    # Write before evaluating either held-out energy.
    output.with_name(output.stem+'-predictions.json').write_text(json.dumps(record,indent=2)+'\n')
    for p in predictions:
        e,t=p['energy'],p['predicted_t']
        p['below_C']=calc(e,.9*t);p['above_C']=calc(e,1.1*t)
        p['actual_t']=brentq(lambda x:calc(e,x),.06,.25,xtol=1e-5)
        p['relative_error']=(t-p['actual_t'])/p['actual_t']
        print('Check',p,flush=True)
    record['rows']=rows
    output.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--grid',nargs=3,type=int,default=[21,20,32])
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();run(a.out,a.grid)
