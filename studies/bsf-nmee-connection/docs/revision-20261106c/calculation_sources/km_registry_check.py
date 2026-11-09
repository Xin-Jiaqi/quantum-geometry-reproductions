from runtime import STUDY, OUTPUT
from pathlib import Path
import sys,csv,numpy as np
r=STUDY;sys.path[:0]=[str(r/'scripts'),str(r),str(r.parents[1]/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from run_fixed_stack_layer_cxx import response, fixed_stack_matrices
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
out=OUTPUT
k,w,_=integration_grid(AuditBilayer(),41,45,72)
with (out/'registry_robustness.csv').open('w') as f:
 writer=None
 for lam,dy in [(x,0) for x in [.15,.20,.25,.30,.35]]+[(.25,-.02),(.25,.02)]:
  u=np.array(ab_ba_shift(lam))+[dy/np.sqrt(3),dy/np.sqrt(3)]
  # The explicit additional displacement is documented in 60-degree fractional coordinates.
  model=AuditBilayer(AuditSpec(shift=tuple(u),staggered_eV=.045,tunneling_eV=1.))
  p=[fixed_stack_matrices(model,k[i:i+512]) for i in range(0,len(k),512)];h=np.concatenate([v[0] for v in p]);v=np.concatenate([v[1] for v in p]);del p
  for t,U in [(10,0),(20,0),(10,12)]:
   row={'lambda':lam,'extra_fraction':dy,'shift60_1':u[0],'shift60_2':u[1],'mass_meV':45,'grid':'41x45x72','step_T':.0625};row.update(response(model,h,v,w,t,U,.0625))
   if writer is None:writer=csv.DictWriter(f,fieldnames=row,lineterminator='\n');writer.writeheader()
   writer.writerow(row);f.flush();print(lam,dy,t,U,row['bottom'],row['top'],row['total'],row['gap_meV'],flush=True)
