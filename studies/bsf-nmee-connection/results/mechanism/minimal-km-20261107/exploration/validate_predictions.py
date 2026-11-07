from pathlib import Path
import json,sys
import numpy as np
repo=Path('/Users/jiaqi/Desktop/worksofJiaqi/quantum-geometry-reproductions')
sys.path.insert(0,str(repo/'studies/bsf-nmee-connection/scripts'))
from test_remote_empty_mechanism import AuditBilayer,integration_grid,ingredients,selective_extra,response
p=Path(__file__).parent;pred=json.loads((p/'main-thresholds-predictions.json').read_text())
k,w,_=integration_grid(AuditBilayer(),81,80,120);rows=[]
for r in pred['predictions']:
 for factor in [.9,1.1]:
  c=np.zeros(3);vmax=-np.inf;cmin=np.inf;t=r['predicted_t']*factor
  for start in range(0,len(k),512):
   data=ingredients(k[start:start+512]);h,v=selective_extra(data,r['energy'],.01,t,.1)
   val,ev,ec=response(h,v,data[-3],data[-2],w[start:start+512],g=data[-1]);c+=val
   vmax=max(vmax,ev);cmin=min(cmin,ec)
  rows.append(dict(energy=r['energy'],t=t,factor=factor,C=c.tolist(),gap=cmin-vmax))
  print(rows[-1],flush=True)
(p/'heldout-fine.json').write_text(json.dumps(dict(grid=[81,80,120],rows=rows),indent=2)+'\n')
