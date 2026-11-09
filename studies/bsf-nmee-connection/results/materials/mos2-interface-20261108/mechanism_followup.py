"""Distinguish mixing of electronic states from gap changes, using existing references."""
from pathlib import Path
import os,socket,json,numpy as np
from analyze import Response,embedding,adj
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
d=np.load('A/cache.npz');out={}
for nb in (96,112):
 r=Response(d,nb)
 for frozen in (4,8,12):
  lift,_=embedding(d,frozen,nb);h=adj(lift)@r.h@lift;lower=d['orb_bottom'];inter=lower[:,None]!=lower[None,:]
  t=h*inter;ht=lift@t@adj(lift);weak=r.h-ht;ew,uw=np.linalg.eigh(weak);es,us=np.linalg.eigh(r.h)
  rec={'spectral':{'weak_states_strong_energies':r.response((uw*es[:,None,:])@adj(uw)),
                   'strong_states_weak_energies':r.response((us*ew[:,None,:])@adj(us))},'channels':{}}
  pv=np.zeros_like(h)
  for mask in (lower,~lower):
   ids=np.where(mask)[0];el,ul=np.linalg.eigh(h[:,ids][:,:,ids]);v=ul[:,:,:14]
   block=v@adj(v)
   for i,ii in enumerate(ids):pv[:,ii,ids]=block[:,i,:]
  pc=np.eye(h.shape[-1])-pv
  vv=pv@t@pv;cc=pc@t@pc;vc=pv@t@pc+pc@t@pv
  assert abs(vv+cc+vc-t).max()<1e-8
  for name,term in [('valence_valence',vv),('conduction_conduction',cc),('valence_conduction',vc),('same_sector',vv+cc)]:
   rec['channels'][name]=r.response(weak+lift@term@adj(lift))
  out[f'{nb}_{frozen}']=rec
  Path('mechanism_followup.json').write_text(json.dumps(out,indent=2))
  print(nb,frozen,rec,flush=True)
