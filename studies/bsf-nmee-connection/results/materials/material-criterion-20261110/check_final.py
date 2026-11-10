from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'summary.json').read_text());ref=json.loads((ROOT/'refine.json').read_text());assert ref['complete']
assert ref['inputs']==s['inputs'];assert all(r['valid'] for r in ref['rows'])
ranges=[]
keys=sorted({(r['material'],r['nb'],r['frozen_empty'],r['mode'],r['alpha']) for r in ref['rows']})
for mat,nb,f,mode,a in keys:
 rr={r['case']:r for r in ref['rows'] if (r['material'],r['nb'],r['frozen_empty'],r['mode'],r['alpha'])==(mat,nb,f,mode,a)}
 g=(np.array(rr['plus']['C_total_lower_upper'])-rr['minus']['C_total_lower_upper'])/(rr['plus']['control']-rr['minus']['control'])
 old=next(r for r in s['rows'] if (r['material'],r['nb'],r['frozen_empty'],r['mode'],r['alpha'])==(mat,nb,f,mode,a))
 gg=np.array([old['G_'+l] for l in ['net','lower','upper']]);ranges.append({'material':mat,'mode':mode,'alpha':a,'G_refined':g.tolist(),'G_previous':gg.tolist(),'abs_error':abs(g-gg).tolist(),'relative_error':(abs(g-gg)/np.maximum(abs(g),1e-12)).tolist(),'all_signs_preserved':bool(np.all(g*gg>0))})
assert all(r['all_signs_preserved'] for r in ranges)
# Confirm physical full screened responses against the independent preceding kernel.
archive=ROOT.parent
oldpath=archive/'material-gate-split-20261110/result.json'
old=json.loads(oldpath.read_text());errs=[]
for r in s['rows']:
 if r['alpha']==1 and r['mode']=='screened_external':
  q=next(x for x in old['rows'] if x['material']==r['material'] and x['nb']==r['nb'] and x['b_eV']==.002 and x['matching']=='anchor')
  errs.extend(abs(r['G_'+l]-q['total'][i]) for i,l in enumerate(['net','lower','upper']))
assert max(errs)<1e-6
out={'refinement':ranges,'all_signs_preserved':True,'max_layer_relative_change':max(x for r in ranges for x in r['relative_error'][1:]),'physical_screened_G_vs_prior_max_error':max(errs),'min_gap_eV':min(r['gap_eV'] for r in ref['rows']),'mixing_diagnostics':ref['mixing_diagnostics'],'passed':True}
(ROOT/'final_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
