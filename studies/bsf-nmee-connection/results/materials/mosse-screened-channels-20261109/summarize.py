"""Tabulate the joint response/intervention test without selecting outcomes."""
from pathlib import Path
import json,csv
import numpy as np
root=Path(__file__).resolve().parent
d=json.loads((root/'result.json').read_text());assert d['complete']
groups={}
for row in d['rows']:
    key=(row['nb'],row['frozen_empty'],row['channel'],row['alpha'])
    groups.setdefault(key,{})[row['case']]=row
table=[]
for (nb,frozen,channel,alpha),rows in groups.items():
    assert set(rows)=={'zero','fm010','fp010'}
    c=np.array(rows['zero']['C_total_lower_upper'])
    g=(np.array(rows['fp010']['C_total_lower_upper'])-rows['fm010']['C_total_lower_upper'])/.2
    r={'nb':nb,'frozen_empty':frozen,'channel':channel,'alpha':alpha}
    r.update({f'C0_{label}':float(x) for label,x in zip(('net','lower','upper'),c)})
    r.update({f'G_{label}':float(x) for label,x in zip(('net','lower','upper'),g)})
    r['upper_absolute_share']=float(abs(g[2])/(abs(g[1])+abs(g[2])))
    r['min_gap_eV']=min(v['gap_eV'] for v in rows.values());table.append(r)
with (root/'summary.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(table[0]),lineterminator='\n');writer.writeheader();writer.writerows(table)
full_by_case={}
for case in ('fm010','fp010'):
    old=json.loads((root.parent/'mosse-gate-scf-20261109'/case/'result.json').read_text())
    full_by_case[case]=max(float(np.max(abs(np.array(r['C_total_lower_upper'])-
        old['materials']['A']['response'][str(r['nb'])]['standard']['C_total_lower_upper'])))
        for r in d['rows'] if r['case']==case and r['channel']=='full')
    assert full_by_case[case]<1e-8, 'Original SCF baseline not reproduced'
checks={'complete':d['complete'],'rows':len(d['rows']),
    'baseline_max_errors':full_by_case,
    'min_gap_eV':min(r['gap_eV'] for r in d['rows']),
    'max_transport_deviation':max(r['representation_max_singular_deviation'] for r in d['checks']),
    'max_half_probe_change':max(r.get('half_probe_max_change',0) for r in d['rows']),
    'max_even_chi':max(max(abs(x) for x in r['even_chi']) for r in d['rows'])}
(root/'audit.json').write_text(json.dumps(checks,indent=2));print(checks)
for r in table:
    if r['nb']==112 and r['frozen_empty']==8:print(r)
