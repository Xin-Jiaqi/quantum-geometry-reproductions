"""Summarize completed calculations only; never fabricate pending results."""
import csv,json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parent
rows=[]
for p in sorted((root/'results').glob('*.json')):
 d=json.loads(p.read_text());c=np.array(d['C_xyz_total_bottom_top_1e-9']);err=d['B_even_residual_xyz']
 rows.append({'structure':d['structure'],'gate_Vnm':d['gate_physical_Vnm'],'mesh':d['mesh'],'probe_VA':d['probe_VA'],'zeeman_eV':d['zeeman_eV'],'C_total':c[0,0],'C_bottom':c[0,1] if d['structure']!='mono' else None,'C_top':c[0,2] if d['structure']!='mono' else None,'B_even_total':None if err is None else err[0][0],'gap_min_eV':min(d['gaps_eV'].values()),'source':p.name})
if not rows:raise SystemExit('No completed material responses yet')
with (root/'summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
lookup={(d['structure'],d['gate_Vnm']):d for d in rows if d['mesh']==6 and d['probe_VA']==.004 and d['zeeman_eV']==.016}
checks={}
if ('A',0) in lookup and ('B',0) in lookup:
 a,b=lookup['A',0],lookup['B',0]
 checks['mirror_sum_total']=a['C_total']+b['C_total'];checks['mirror_sum_bottomA_topB']=a['C_bottom']+b['C_top'];checks['mirror_sum_topA_bottomB']=a['C_top']+b['C_bottom']
if ('mono',.25) in lookup and ('mono',-.25) in lookup:checks['mono_gate_oddness_residual']=lookup['mono',.25]['C_total']+lookup['mono',-.25]['C_total']
(root/'summary.json').write_text(json.dumps({'rows':rows,'checks':checks,'units':'1e-9 (nm V T)^-1, Pauli spin density','missing_rows_are_pending_not_zero':True},indent=2))
print(json.dumps({'rows':rows,'checks':checks},indent=2))
