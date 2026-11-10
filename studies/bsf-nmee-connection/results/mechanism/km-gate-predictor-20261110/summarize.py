"""Audit zero-bias predictions against held-out physical offsets."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(name):
 d=json.loads((ROOT/name).read_text());assert d['complete'];return d
pred=read('prediction-fine.json');obs=read('verification-fine.json')
mirror=read('mirror-zero.json');mv=read('mirror-verification.json');refined=read('prediction-refined.json')
sealed=json.loads((ROOT/'mirror-prediction.json').read_text())
assert sealed['source_sha256']==hashlib.sha256((ROOT/'prediction-fine.json').read_bytes()).hexdigest()
audit=[];deltas=[]
for r in pred['rows']:
 t=r['t_meV'];c0=np.array(r['chi']);g=np.array(r['tangent']);states={x['U_meV']:x for x in obs['rows'] if x['t_meV']==t}
 for u in (.25,.5):
  fd=(np.array(states[u]['chi'])-states[-u]['chi'])/(2*u/1000)
  audit.append(dict(t_meV=t,half_window_meV=u,G_analytic=g.tolist(),G_finite=fd.tolist(),relative_error=(abs(fd-g)/np.maximum(abs(g),1e-12)).tolist()))
 for u,x in states.items():
  actual=np.array(x['chi'])-c0;expected=g*u/1000
  deltas.append(dict(t_meV=t,U_meV=u,delta_C=actual.tolist(),linear_prediction=expected.tolist(),sign_match=(np.sign(actual)==np.sign(expected)).tolist(),one_sided_relative_error=(abs(actual-expected)/np.maximum(abs(actual),1e-12)).tolist()))
for r in pred['rows']+mirror['rows']+refined['rows']:
 assert abs(sum(r['chi'][:2])-r['chi'][2])<1e-10
 assert abs(sum(r['tangent'][:2])-r['tangent'][2])<1e-8
 assert np.max(abs(np.array(r['tangent'])-r['energy_shift']-np.array(r['matrix_change'])))<1e-8
mr=mirror['rows'][0]
mirror_error=max(np.max(abs(np.array(mr['chi'])-sealed['C_zero'])),np.max(abs(np.array(mr['tangent'])-sealed['G_zero'])))
a0=pred['rows'][0];b0=mr
local_change={}
for label,zero,rows in [('A',a0,[x for x in obs['rows'] if x['t_meV']==10]),('B',b0,mv['rows'])]:
 local_change[label]={str(r['U_meV']):(abs(np.array(r['chi'][:2]))-abs(np.array(zero['chi'][:2]))).tolist() for r in rows if r['U_meV']>0}
summary={'order':['lower','upper','net'],'central_derivative_checks':audit,'finite_bias_predictions':deltas,'all_direction_predictions_correct':all(all(r['sign_match']) for r in deltas),
 'max_central_relative_error':max(max(x['relative_error']) for x in audit),'max_one_sided_relative_error':max(max(x['one_sided_relative_error']) for x in deltas),
 'mirror_absolute_error':float(mirror_error),'mirror_positive_offset_magnitude_changes':local_change,
 'minimum_common_gap_meV':min(r['common_gap_meV'] for d in [pred,obs,mirror,mv,refined] for r in d['rows']),
 'refined_t10_upper_share':abs(refined['rows'][0]['tangent'][1])/sum(abs(x) for x in refined['rows'][0]['tangent'][:2]),
 'limits':['Full zero-bias response derivative, not an independently established few-parameter universal rule.',
 'Energy-only net sign is mesh-sensitive in coarse data; no robust wrong-net-sign claim.',
 'Small lower slope has large one-sided nonlinear relative errors even within selected offsets.',
 'Finite-B eigenbasis diagnostics assume resolved bands; no band-by-band interpretation at exact degeneracy.',
 'KM internal U excludes intralayer field variation and self-consistent screening.',
 'Material fixed-subspace channel interventions are not these instantaneous eigenstate-rotation sectors.']}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
rows=[]
for label,d in [('prediction',pred),('verification',obs),('mirror_zero',mirror),('mirror_verification',mv),('refined',refined)]:
 for r in d['rows']:
  for i,layer in enumerate(('lower','upper','net')):
   rows.append(dict(dataset=label,fraction=r['fraction'],t_meV=r['t_meV'],U_meV=r['U_meV'],layer=layer,C=r['chi'][i],G=r.get('tangent',[None]*3)[i],G_energy=r.get('energy_shift',[None]*3)[i],G_matrix=r.get('matrix_change',[None]*3)[i],gap_meV=r['common_gap_meV']))
with (ROOT/'source_data.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps({k:v for k,v in summary.items() if k not in ['central_derivative_checks','finite_bias_predictions']},indent=2))
