import sys,json
from pathlib import Path
import numpy as np
repo=Path('/Users/jiaqi/Desktop/worksofJiaqi/quantum-geometry-reproductions')
sys.path[:0]=[str(repo/'studies/bsf-nmee-connection'),str(repo/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.remote_empty_km import ingredients,assemble,response
from bsf_nmee.bilayer_audit import AuditBilayer,integration_grid
k,w,_=integration_grid(AuditBilayer(),21,20,32)
rows=[]
for site in [0,1]:
 data=ingredients(k,site=site)
 for energy in [.06,.12,.25,.5]:
  for t in [0.,.01,.03,.06,.12]:
   for channel in ['all','VC']:
    h,v=assemble(data,energy,0,t,channel)
    c,ev,ec=response(h,v,data[-3],data[-2],w,g=data[-1])
    r=dict(site=site,energy=energy,t=t,channel=channel,C=c.tolist(),gap=ec-ev)
    rows.append(r)
    if t==0 or c[0]>0:print(r,flush=True)
 Path(__file__).with_name('explore.json').write_text(json.dumps(rows,indent=2)+'\n')
print('DONE',len(rows),flush=True)
