"""Full-space selected interlayer interventions; never additive contributions.
Run once each independently strained DFT cache is available, on PBS only.
"""
from pathlib import Path
import os,socket,json
import numpy as np
from analyze import Response,embedding,adj
from mixing_channels import split_channels
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();label=os.environ['CASE'];d=np.load(root/label/'A/cache.npz')
bottom=d['orb_bottom'];layers=[np.flatnonzero(bottom),np.flatnonzero(~bottom)]
assert len(layers[0])==len(layers[1])==22
result={'case':label,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'units':'1e-9 (nm V T)^-1','scope':'Fixed zero-field layer projectors, 14 occupied of 22 d+p states per layer; interventions retain all KS response vertices and unselected states. Not additive electrostatic/hybridization decomposition.','runs':{}}
combinations={'none':(), 'VV_CC':('VV','CC'),'VC':('VC',),'all':('VV','CC','VC')}
for nb in (96,112):
 r=Response(d,nb)
 for frozen in (4,8,12):
  lift,s=embedding(d,frozen,nb);hloc=adj(lift)@r.h@lift
  intra,parts,p,checks=split_channels(hloc,layers,14)
  embedded={k:lift@v@adj(lift) for k,v in parts.items()};weak=r.h-sum(embedded.values())
  rows={}
  for name,include in combinations.items():
   h=weak+sum((embedded[k] for k in include),np.zeros_like(weak))
   try: rows[name]=r.response(h)
   except ValueError as err: rows[name]={'invalid_control':str(err)}
  # Keep the zero-field valence/empty subspaces fixed in all probe calculations.
  norms={k:float(np.sqrt(np.mean(np.sum(abs(v)**2,axis=(1,2))))) for k,v in parts.items()}
  record={'min_projection_singular':s,'checks':checks,'rows':rows,'selected_channel_norm_eV':norms}
  if frozen==8:
   record['mixing_scan']=[dict(alpha=a,**r.response(r.h+(a-1)*sum(embedded.values()))) for a in (.5,.75,1)]
  result['runs'][f'{nb}_{frozen}']=record
  (root/label/'mechanism.json').write_text(json.dumps(result,indent=2));print(label,nb,frozen,{k:round(v['C_total_lower_upper'][1],5) if 'C_total_lower_upper' in v else v for k,v in rows.items()},flush=True)
