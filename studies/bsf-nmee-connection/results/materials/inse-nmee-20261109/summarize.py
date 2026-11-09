"""Physical layers, paired gates and mirror checks, retaining every band cutoff."""
from pathlib import Path
import json,csv
import numpy as np
p=Path(__file__).resolve().parent;names=['mono0','monom','monop','A0','B0','Am','Ap']
data={n:json.loads((p/n/'result.json').read_text()) for n in names if (p/n/'result.json').exists()}
rows=[]
for name,d in data.items():
    for nb,ans in d['response'].items():
        v=ans['standard'];c=v['C_total_lower_upper'];half=ans['half_b']['C_total_lower_upper']
        row={'name':name,'nb':int(nb),'F':d['F_external_V_nm'],'C_net':c[0],
            'C_lower':None if d['label']=='mono' else c[1],
            'C_upper':None if d['label']=='mono' else c[2],
            'gap_eV':v['gap_eV'],'probe_half_max_physical_change':abs(c[0]-half[0]) if d['label']=='mono' else float(max(abs(np.array(c)-half))),
            'even_chi_net':v['even_chi'][0]}
        rows.append(row)
if rows:
    with (p/'points.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
slopes=[]
for stem,plus,minus in [('A','Ap','Am'),('mono','monop','monom')]:
    if plus not in data or minus not in data:continue
    for nb in data[plus]['response']:
        cp=np.array(data[plus]['response'][nb]['standard']['C_total_lower_upper'])
        cm=np.array(data[minus]['response'][nb]['standard']['C_total_lower_upper']);g=(cp-cm)/.2
        slopes.append({'system':stem,'nb':int(nb),'G_net':float(g[0]),
            'G_lower':None if stem=='mono' else float(g[1]),'G_upper':None if stem=='mono' else float(g[2]),
            'matched_pair_G':float(2*g[0]) if stem=='mono' else None})
if slopes:
    with (p/'slopes.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(slopes[0]),lineterminator='\n');w.writeheader();w.writerows(slopes)
checks={'complete_points':list(data),'all_seven_complete':len(data)==7}
if 'mono0' in data:checks['mono_zero_C']={nb:v['standard']['C_total_lower_upper'][0] for nb,v in data['mono0']['response'].items()}
if 'monop' in data and 'monom' in data:checks['mono_gate_odd_residual']={nb:data['monop']['response'][nb]['standard']['C_total_lower_upper'][0]+data['monom']['response'][nb]['standard']['C_total_lower_upper'][0] for nb in data['monop']['response']}
if 'A0' in data and 'B0' in data:
    checks['mirror_residuals']={}
    for nb in data['A0']['response']:
        a=np.array(data['A0']['response'][nb]['standard']['C_total_lower_upper'])
        b=np.array(data['B0']['response'][nb]['standard']['C_total_lower_upper'])
        checks['mirror_residuals'][nb]=(a+b[[0,2,1]]).tolist()
if rows:checks['minimum_gap_eV']=min(r['gap_eV'] for r in rows)
(p/'summary.json').write_text(json.dumps({'checks':checks,'points':rows,'slopes':slopes},indent=2))
print(json.dumps(checks,indent=2))
for r in rows:print(r)
for r in slopes:print(r)
