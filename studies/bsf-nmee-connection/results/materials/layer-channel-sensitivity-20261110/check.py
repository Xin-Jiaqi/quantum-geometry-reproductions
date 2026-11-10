"""Cross-check new summary against independent previously archived totals."""
from pathlib import Path
import json,argparse,hashlib
import numpy as np
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--materials-root',type=Path,default=root.parent);args=p.parse_args()
s=json.loads((root/'summary.json').read_text());prior=json.loads((args.materials_root/'material-gate-split-20261110/result.json').read_text())
errors=[]
for r in s['rows']:
 if r['channel']!='full':continue
 old=next(x for x in prior['rows'] if x['material']==r['material'] and x['nb']==r['nb'] and x['b_eV']==.002 and x['matching']=='anchor')
 for i,l in enumerate(('net','lower','upper')):
  errors.extend([abs(r[f'G_{l}']-old['total'][i]),abs(r[f'Cminus_{l}']-old['Cminus'][i]),abs(r[f'Cplus_{l}']-old['Cplus'][i])])
assert max(errors)<1e-6
for source in s['sources']:
 assert hashlib.sha256((args.materials_root/source['path']).read_bytes()).hexdigest()==source['sha256']
ranges={(r['material'],r['channel'],r['alpha']):r for r in s['ranges']}
a=ranges['MoSSe','Vb_Ct',0];b=ranges['InSe','Vb_Ct',0]
assert a['all_lower_zero_gate_sign_changes'] and not b['all_lower_zero_gate_sign_changes']
assert a['ranges']['R_ratio_to_full_upper'][1]<1 and b['ranges']['R_ratio_to_full_upper'][0]>1
assert not a['any_endpoint_layer_sign_change'] and not b['any_endpoint_layer_sign_change']
for material,n in [('MoSSe',6),('InSe',9)]:
 assert s['checks'][material]['row_count']==7*n
assert len(s['rows'])==105
out={'full_response_vs_previous_independent_kernel_max_error':max(errors),'source_hashes_verified':True,'all_105_groups_retained':True,'primary_directional_channel_claim_holds_all_windows':True,'primary_controls_no_gate_endpoint_layer_sign_change':True,'passed':True}
(root/'checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
