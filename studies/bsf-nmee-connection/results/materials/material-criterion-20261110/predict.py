"""Freeze alpha-Taylor predictions before held-out jobs; never refit on validation."""
from pathlib import Path
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parent
raw=ROOT/'train.json';d=json.loads(raw.read_text());assert d['complete']
groups=sorted({(r['material'],r['nb'],r['frozen_empty'],r['mode']) for r in d['rows']})
out={'saved_before_validation':True,'training_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'training_alpha':[-.25,0,.25],
 'held_out_new_alpha':[-.5,.5],'physical_alpha_1':'Previously known in some modes; retrospective transfer/extrapolation check, not blind.',
 'method':'Local alpha polynomial around removed selected channel; not an assumed weak-total-mixing expansion.', 'rows':[],'invalid_groups':[]}
for mat,nb,f,mode in groups:
 rs=[r for r in d['rows'] if (r['material'],r['nb'],r['frozen_empty'],r['mode'])==(mat,nb,f,mode)]
 if not all(r['valid'] for r in rs):out['invalid_groups'].append([mat,nb,f,mode]);continue
 def get(a,case):return next(r for r in rs if r['alpha']==a and r['case']==case)
 vals={}
 for a in [-.25,0,.25]:
  minus,plus=get(a,'minus'),get(a,'plus');delta=plus['control']-minus['control']
  vals[a]=(np.array(plus['C_total_lower_upper'])-minus['C_total_lower_upper'])/delta
 A=vals[0];B=(vals[.25]-vals[-.25])/.5;D=(vals[.25]+vals[-.25]-2*A)/(2*.25**2)
 out['rows'].append({'material':mat,'nb':nb,'frozen_empty':f,'mode':mode,'A':A.tolist(),'B':B.tolist(),'D':D.tolist(),
 'predictions':[{'alpha':a,'G':(A+B*a+D*a*a).tolist(),'abs_effect_vs_removed':(abs(A+B*a+D*a*a)-abs(A)).tolist()} for a in [-.5,.5,1.]]})
(ROOT/'predictions.json').write_text(json.dumps(out,indent=2)+'\n');print('Frozen predictions:',len(out['rows']),'invalid',len(out['invalid_groups']))
