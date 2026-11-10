"""Evaluate fixed predictions; failures and invalid controls remain in output."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
train=json.loads((ROOT/'train.json').read_text());val=json.loads((ROOT/'validate.json').read_text());pred=json.loads((ROOT/'predictions.json').read_text())
assert train['complete'] and val['complete'];assert pred['training_sha256']==hashlib.sha256((ROOT/'train.json').read_bytes()).hexdigest();assert train['inputs']==val['inputs']
rows=[];invalid=[]
for source,d in [('train',train),('validate',val)]:
 keys=sorted({(r['material'],r['nb'],r['frozen_empty'],r['mode'],r['alpha']) for r in d['rows']})
 for mat,nb,f,mode,a in keys:
  rs={r['case']:r for r in d['rows'] if (r['material'],r['nb'],r['frozen_empty'],r['mode'],r['alpha'])==(mat,nb,f,mode,a)}
  if len(rs)!=3 or not all(r['valid'] for r in rs.values()):invalid.append({'source':source,'key':[mat,nb,f,mode,a],'states':rs});continue
  cm=np.array(rs['minus']['C_total_lower_upper']);cp=np.array(rs['plus']['C_total_lower_upper']);c0=np.array(rs['zero']['C_total_lower_upper']);g=(cp-cm)/(rs['plus']['control']-rs['minus']['control'])
  r={'source':source,'material':mat,'nb':nb,'frozen_empty':f,'mode':mode,'alpha':a,'min_gap_eV':min(x['gap_eV'] for x in rs.values())}
  for i,l in enumerate(['net','lower','upper']):r.update({f'C0_{l}':c0[i],f'G_{l}':g[i],f'gate_endpoint_sign_change_{l}':bool(cm[i]*cp[i]<=0)})
  rows.append(r)
tests=[]
for r in rows:
 if r['source']!='validate':continue
 key=(r['material'],r['nb'],r['frozen_empty'],r['mode']);pr=next((x for x in pred['rows'] if (x['material'],x['nb'],x['frozen_empty'],x['mode'])==key),None)
 if pr is None:continue
 ex=next(x for x in pr['predictions'] if x['alpha']==r['alpha']);actual=np.array([r['G_'+l] for l in ['net','lower','upper']]);expected=np.array(ex['G']);baseline=np.array(pr['A'])
 test={k:r[k] for k in ['material','nb','frozen_empty','mode','alpha']};test['blind_new_control']=r['alpha']!=1
 for i,l in enumerate(['net','lower','upper']):test.update({f'actual_{l}':actual[i],f'predicted_{l}':expected[i],f'baseline_{l}':baseline[i],f'sign_correct_{l}':bool(actual[i]*expected[i]>0),f'absolute_error_{l}':abs(actual[i]-expected[i]),f'relative_error_{l}':abs(actual[i]-expected[i])/max(abs(actual[i]),1e-12),f'gain_direction_correct_{l}':bool((abs(expected[i])-abs(baseline[i]))*(abs(actual[i])-abs(baseline[i]))>0)})
 tests.append(test)
for file,records in [('all_controls.csv',rows),('prediction_tests.csv',tests)]:
 if records:
  with (ROOT/file).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)
out={'rows':rows,'tests':tests,'invalid':invalid,'inputs':train['inputs'],'claim':'Quadratic transfer is a hypothesis tested without post-validation refitting. No universal criterion is assumed.','summary':[]}
for mat in ['MoSSe','InSe','MoS2']:
 for mode in ['internal_step','screened_external']:
  rr=[x for x in tests if x['material']==mat and x['mode']==mode]
  if not rr:continue
  for kind in ['new_controls','physical_1']:
   subset=[x for x in rr if (x['alpha']!=1)==(kind=='new_controls')]
   out['summary'].append({'material':mat,'mode':mode,'kind':kind,'n_cases':len(subset),**{f'{l}_sign_correct':sum(x[f'sign_correct_{l}'] for x in subset) for l in ['lower','upper','net']},'max_layer_relative_error':max(x[f'relative_error_{l}'] for x in subset for l in ['lower','upper'])})
(ROOT/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['summary'],indent=2));print('invalid',len(invalid))
