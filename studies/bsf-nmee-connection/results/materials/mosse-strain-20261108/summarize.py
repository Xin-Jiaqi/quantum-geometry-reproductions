"""Local summary of completed results. Does not infer missing strain points."""
from pathlib import Path
import json,csv
root=Path(__file__).resolve().parent
labels=('em02','e000','ep02');rows=[];checks=[];mechanism=[]
for label in labels:
 base=root/label;d=json.loads((base/'result.json').read_text());assert d['stage']=='complete'
 strain=d['strain'];geom=d['geometry']
 assert geom['max_z_force_eVA']<=.012 and geom['max_inplane_force_eVA']<=.012
 for material in ('A','mono'):
  exported=json.loads((base/material/'export.json').read_text())
  assert exported['metric_max']<1e-4 and exported['wave_link_gauge_error']<.08
  for nb,rec in d[material+'_response'].items():
   standard=rec['standard'];half=rec['half_b'];assert standard['gap_eV']>.05
   c=standard['C_total_lower_upper'];ch=half['C_total_lower_upper']
   rows.append({'strain_pct':100*strain,'system':material,'nb':nb,'reference_end_A':'','lower':c[1] if material=='A' else c[0],'upper':c[2] if material=='A' else '', 'total':c[0],'gap_eV':standard['gap_eV']})
   checks.append({'case':label,'system':material,'nb':nb,'max_b_step_change':max(abs(x-y) for x,y in zip(c,ch)),'min_neighbor_singular':standard['min_neighbor_singular'],'SS_gap_A':geom['interface_SS_gap_A']})
 ref=json.loads((base/'reference.json').read_text())
 for r in ref['rows']:
  rows.append({'strain_pct':100*strain,'system':'electrostatic','nb':r['nb'],'reference_end_A':r['taper_end_A'],'lower':r['C_lower'],'upper':'','total':'','gap_eV':r['gap_eV']})
 m=json.loads((base/'mechanism.json').read_text());assert len(m['runs'])==6
 for key,r in m['runs'].items():
  for name,v in r['rows'].items():
   mechanism.append({'strain_pct':100*strain,'window':key,'control':name,'lower':v.get('C_total_lower_upper',[None,None])[1],'gap_eV':v.get('gap_eV'),'invalid':v.get('invalid_control','')})
for name,data in [('responses.csv',rows),('mechanism.csv',mechanism)]:
 with (root/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
(root/'checks.json').write_text(json.dumps(checks,indent=2))
for strain in (-2.,0.,2.):
 rr=[r for r in rows if r['strain_pct']==strain]
 print(strain,{kind:(min(r['lower'] for r in rr if r['system']==kind),max(r['lower'] for r in rr if r['system']==kind)) for kind in ('A','mono','electrostatic')})
