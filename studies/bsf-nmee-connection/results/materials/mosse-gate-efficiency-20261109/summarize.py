"""Summarize all computed gate curves without normalizing by a small C(0)."""
from pathlib import Path
import json,csv
from collections import defaultdict
import numpy as np
root=Path(__file__).resolve().parent
d=json.loads((root/'scan.json').read_text());assert d.get('complete')
rows=[];groups=defaultdict(dict);invalid=[]
for source,records in [('bilayer',d['bilayer']),('reference',d['references'])]:
    for r in records:
        key=(r.get('control',r.get('kind')),r['nb'],r.get('frozen_empty',''),r.get('taper_end_A',''))
        if 'invalid' in r:
            invalid.append({'key':key,'field':r['F_model'],'detail':r['invalid']});continue
        c=r.get('C',r.get('C_total_lower_upper'))
        assert abs(c[0]-c[1]-c[2])<1e-7
        groups[key][r['F_model']]=r
        rows.append(dict(kind=key[0],nb=key[1],frozen_empty=key[2],taper_end_A=key[3],
                         F_model=r['F_model'],C_net=c[0],C_lower=c[1],C_upper=c[2],gap_eV=r['gap_eV']))
slopes=[]
for key,g in groups.items():
    if 0. not in g:continue
    def c(r):return np.array(r.get('C',r.get('C_total_lower_upper')))
    for delta in (.025,.05,.1,.25,.5):
        if delta not in g or -delta not in g:continue
        gg=(c(g[delta])-c(g[-delta]))/(2*delta)
        slopes.append(dict(kind=key[0],nb=key[1],frozen_empty=key[2],taper_end_A=key[3],
                           half_interval=delta,G_net=gg[0],G_lower=gg[1],G_upper=gg[2]))
for name,data in [('curves.csv',rows),('slopes.csv',slopes)]:
    with (root/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
summary={'accepted_curve_rows':len(rows),'invalid_rows':invalid,'central_slopes':{},
         'min_gap_eV':min(r['gap_eV'] for r in rows),
         'max_half_B_change':max(r.get('half_probe_step_max_change',0) for r in d['bilayer']),
         'field_calibration':'imposed one-electron model slope; no external screening calibration'}
for kind in sorted(set(s['kind'] for s in slopes)):
    rr=[s for s in slopes if s['kind']==kind and s['half_interval']==.025]
    summary['central_slopes'][kind]={a:[min(s[a] for s in rr),max(s[a] for s in rr)] for a in ['G_net','G_lower','G_upper']}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
