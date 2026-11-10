from pathlib import Path
import json,sys
import numpy as np
import run as km
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'summary.json').read_text())
k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
plus=km.AuditBilayer(km.AuditSpec(shift=km.ab_ba_shift(.25),staggered_eV=.045,tunneling_eV=.002))
minus=km.AuditBilayer(km.AuditSpec(shift=km.ab_ba_shift(.25),staggered_eV=.045,tunneling_eV=-.002))
l=km.P_TOP-km.P_BOTTOM
gauge=max(float(abs(minus.matrices(k,axis=a)-l@plus.matrices(k,axis=a)@l).max()) for a in [None,0,1])
assert gauge<1e-12
assert all(float(abs(l@sp@l-sp).max())<1e-12 for sp in plus.spins.values())
assert all(e['signs_correct'] for e in s['weak_validation']) and s['phase_signs_correct']
assert s['strong_contrast'][0]['full_to_removed_abs_G_ratio']>1 and s['strong_contrast'][1]['full_to_removed_abs_G_ratio']<1
assert s['all_insulating'] and s['layer_sum_G_max']<1e-8
oldroot=ROOT.parent/'km-gate-predictor-20261110'
old=json.loads((oldroot/'prediction-fine.json').read_text())['rows'];comparisons=[]
for r in s['rows']:
 if r['source']=='refined.json' and r['channel']=='full' and r['U_meV']==0:
  q=next(x for x in old if x['t_meV']==r['t_meV'])
  comparisons.append({'t_meV':r['t_meV'],'C_error':max(abs(r['C_'+v]-q['chi'][i]) for i,v in enumerate(['lower','upper','net'])),'G_vs_analytic_error':max(abs(r['G_'+v]-q['tangent'][i]) for i,v in enumerate(['lower','upper','net']))})
assert max(r['G_vs_analytic_error'] for r in comparisons)<1e-5
out={'global_interlayer_sign_gauge_error':gauge,'channel_derivative_audit':km.audit(),'old_kernel_baseline':comparisons,'all_weak_prediction_signs_correct':True,'phase_prediction_signs_correct':True,'new_states':0,'passed':True}
(ROOT/'checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
