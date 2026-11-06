from runtime import STUDY, OUTPUT
from pathlib import Path
import sys,csv,json,numpy as np
from types import SimpleNamespace
r=STUDY;sys.path[:0]=[str(r/'scripts'),str(r),str(r.parents[1]/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from run_fixed_stack_layer_cxx import response
from bsf_nmee.bilayer_audit import AuditBilayer,AuditSpec,integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from bsf_nmee.sliding_tight_binding import P_TOP,P_BOTTOM
out=OUTPUT
mode=sys.argv[1]; grid=(41,45,72) if mode=='robust' else (81,80,120);step=.0625 if mode=='robust' else .03125
cases=[(.24,45),(.25,45),(.26,45)] if mode=='robust' else [(.25,45)]
pts=[(10,U) for U in ([1,2,3,4.5] if mode=='robust' else [1,3,4.5,4.75,5])]
e=np.array([.5,np.sqrt(3)/2]);k,w,_=integration_grid(AuditBilayer(),*grid)
with (out/('gate_extra_'+mode+'.csv')).open('w') as f:
 wr=None
 for lam,mass in cases:
  model=AuditBilayer(AuditSpec(shift=ab_ba_shift(lam),staggered_eV=mass/1000,tunneling_eV=1.))
  h=np.concatenate([model.matrices(k[i:i+512]) for i in range(0,len(k),512)])
  v=np.concatenate([e[0]*model.matrices(k[i:i+512],axis=0)+e[1]*model.matrices(k[i:i+512],axis=1) for i in range(0,len(k),512)])
  spins={name+'_x':e[0]*model.spins[name+'_x']+e[1]*model.spins[name+'_y'] for name in ['bottom','top','total']};spins['total_z']=model.spins['total_z'];probe=SimpleNamespace(mono=model.mono,spins=spins)
  for t,U in pts:
   row={'lambda':lam,'mass_meV':mass,'grid':'x'.join(map(str,grid)),'step_T':step,'probe_angle_deg':60};row.update(response(probe,h,v,w,t,U,step))
   if wr is None:wr=csv.DictWriter(f,fieldnames=row,lineterminator='\n');wr.writeheader()
   wr.writerow(row);f.flush();print(mode,lam,mass,t,U,row['bottom'],row['top'],row['total'],row['common_gap_meV'],flush=True)
