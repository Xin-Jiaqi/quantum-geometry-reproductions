"""Check time-reversal eigenvalue pairing and empty-band sensitivity on PBS."""
import os,socket,json
from pathlib import Path
import numpy as np
if not os.environ.get('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
d=np.load('full-controls-cache.npz'); nk=len(d['k']);n=round(nk**.5)
ks=np.rint(d['k'][:,:2]*n).astype(int)%n
lookup={tuple(k):i for i,k in enumerate(ks)};neg=np.array([lookup[tuple((-k)%n)] for k in ks])
e=np.diagonal(d['diag'],axis1=1,axis2=2).real
out={'node':socket.gethostname(),'job':os.environ['PBS_JOBID'],'energy_TR_max_occupied':float(abs(e[:,:48]-e[neg,:48]).max()),'energy_TR_max_all':float(abs(e-e[neg]).max()),'spin_z_occupied':float(np.trace(d['Z'][:,:48,:48],axis1=1,axis2=2).real.mean()),'spin_x_occupied':float(np.trace(d['spin_total'][:,:48,:48],axis1=1,axis2=2).real.mean()),'cutoffs':[]}
adj=lambda a:a.conj().swapaxes(-1,-2)
for nb in (64,80,96,112):
 h=d['diag'][:,:nb,:nb];z=d['Z'][:,:nb,:nb];links=d['links'][:,:,:nb,:nb]
 spin=d['spin'][:,:nb,:nb];total=d['spin_total'][:,:nb,:nb];P=d['P'][:,:nb,:nb]
 def chi(h,op):
  en,u=np.linalg.eigh(h);occ=u[:,:,:48];emp=u[:,:,48:];deriv=np.zeros_like(occ,dtype=complex)
  for a in range(4):
   v=links[:,a]@occ[d['neighbors'][:,a]]
   deriv+=(-1)**a*n*d['cell'][a//2,0]/(4*np.pi)*v@np.linalg.inv(adj(occ)@v)
  return float((2*((adj(occ)@op@emp)*(adj(emp)@(1j*deriv)).swapaxes(1,2)).real/(en[:,:48,None]-en[:,None,48:])).sum(axis=(1,2)).mean()*6.331927853335898e4)
 rows=[]
 for U in (0.,.15,.6,-.15):
  hu=h+U*(.5*np.eye(nb)-P);r={'U':U}
  for name,op in [('bottom',spin),('total',total)]:
   plus=chi(hu+.002*z,op);minus=chi(hu-.002*z,op)
   r[name]={'C':(plus-minus)/.004,'even':(plus+minus)/2}
  rows.append(r)
 out['cutoffs'].append({'bands':nb,'rows':rows})
Path('zero-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
