"""Check split algebra and same-convention KM limit before material execution."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parent
candidates=[ROOT.parent/'km-gate-predictor-20261110', ROOT.parent.parent/'mechanism/km-gate-predictor-20261110']
km_root=next((p for p in candidates if (p/'predict.py').exists()),None)
if km_root is None: raise FileNotFoundError('Sibling archived KM predictor required')
sys.path.insert(0,str(km_root))
import predict as km
from split_kernel import adj,match,reorder,split
m=km.AuditBilayer(km.AuditSpec(shift=km.ab_ba_shift(.25),staggered_eV=.045,tunneling_eV=.01))
k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
h=m.matrices(k)-.03125*m.mono.g_factor*km.MU_B*m.spins['total_z']
v=.5*m.matrices(k,axis=0)+np.sqrt(3)/2*m.matrices(k,axis=1)
s=.5*m.spins['bottom_x']+np.sqrt(3)/2*m.spins['bottom_y']
e0,u0=np.linalg.eigh(h)
ref=km.kernel(h,v,[s])[0][0]
def fac(e,u):
 sv=(adj(u)@s@u)[:,:4,4:];vv=(adj(u)@v@u)[:,4:,:4].swapaxes(1,2)
 delta=e[:,:4,None]-e[:,None,4:]
 n=-2*np.imag(sv*vv);d=delta**-2
 rr=1j*vv/delta
 assert np.max(abs(n*d-2*np.real(sv*rr)/delta))<1e-10
 return n,d
errors=[]
for du in (1e-5,5e-6):
 es=[]
 for sign in (1,-1):
  e,u=np.linalg.eigh(h+sign*du*km.W)
  perm,_=match(u0,u,np.tile(np.eye(8),(len(k),1,1)),4)
  es.append(fac(*reorder(e,u,perm)))
 ee,mm,tt=split(*es[0],*es[1],2*du)
 errors.append({'du_eV':du,'energy_error':float(np.max(abs(ee.sum((1,2))-ref['energy_shift']))),'matrix_error':float(np.max(abs(mm.sum((1,2))-ref['matrix_change']))),'total_error':float(np.max(abs(tt.sum((1,2))-ref['tangent']))),'closure':float(np.max(abs(tt-ee-mm)))})
assert max(max(v for key,v in x.items() if key!='du_eV') for x in errors)<1e-7
rng=np.random.default_rng(29);perm=np.array([rng.permutation(4).tolist()+(4+rng.permutation(4)).tolist() for _ in k]);ep,up=reorder(e0,u0,perm);up=up*np.exp(1j*rng.random(ep.shape))[:,None,:]
p,_=match(u0,up,np.tile(np.eye(8),(len(k),1,1)),4)
n,d=fac(*reorder(ep,up,p));n0,d0=fac(e0,u0)
phase=float(np.max(abs(n*d-n0*d0)));assert phase<1e-10
out={'same_squared_denominator_KM_limit':errors,'permutation_and_phase_error':phase,'passed':True}
(ROOT/'local_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
