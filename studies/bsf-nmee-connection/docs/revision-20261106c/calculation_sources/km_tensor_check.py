from runtime import STUDY, OUTPUT
from pathlib import Path
import sys,json
from dataclasses import replace
r=STUDY;sys.path[:0]=[str(r),str(r.parents[1]/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift,prepared_B_response,half_filled_spec
from bsf_nmee.sliding_tight_binding import P_TOP,P_BOTTOM
# Apply the already-defined uniform layer bias to H, not to its k derivative.
class Biased(AuditBilayer):
 def __init__(self,spec,U):super().__init__(spec);self.U=U
 def matrices(self,k,axis=None):
  h=super().matrices(k,axis)
  if axis is None:h+=self.U/2000*(P_TOP-P_BOTTOM)
  return h
k,w,_=integration_grid(AuditBilayer(),41,45,72)
out=(OUTPUT/'tensor_check.json');rows=[]
for lam,U in [(0,0),(0,12),(.25,0),(.25,10),(.25,11),(.25,12)]:
 spec=AuditSpec(shift=ab_ba_shift(lam),staggered_eV=.045,tunneling_eV=.010)
 m=Biased(spec,U)
 import numpy as np
 e=np.linalg.eigvalsh(m.matrices(k));lo=e[:,3].max();hi=e[:,4].min();assert hi>lo
 m=Biased(replace(spec,mu_eV=(lo+hi)/2),U)
 row={'lambda':lam,'U_meV':U,'gap_meV':float((hi-lo)*1000),'tensor':{}}
 for obs in ['bottom','top','total']:
  row['tensor'][obs]=prepared_B_response(m,k,w,[0],eta_eV=1e-9,step_T=.0625,observation=obs)['real'][0]
 rows.append(row);out.write_text(json.dumps(rows,indent=2));print(row,flush=True)
