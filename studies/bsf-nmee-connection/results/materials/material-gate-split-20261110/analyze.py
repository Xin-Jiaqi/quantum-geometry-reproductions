"""Summarize every diagnostic, retaining failed mechanism attribution."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
sources=[json.loads((ROOT/f).read_text()) for f in ('result.json','result_refined.json')]
assert all(s['complete'] for s in sources)
assert sources[0]['inputs']==sources[1]['inputs']
rows=sum([s['rows'] for s in sources],[])
fields=['material','nb','b_eV','matching','layer','energy','matrix','total','energy_single_den','matrix_single_den','Cminus','Cplus']
flat=[]
for r in rows:
 for i,layer in enumerate(('net','lower','upper')):
  flat.append({**{k:r[k] for k in fields[:4]},'layer':layer,**{k:r[k][i] for k in fields[5:]}})
with (ROOT/'all_diagnostics.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(flat)
summary={'scope':'Existing 6x6 fixed-ion SCFs, -0.10 to +0.10 V/nm. No new k-mesh or DFT convergence claim.',
 'interpretation':'Stable total changes; energy/matrix partition has no established small-B limit. Cancellation mechanism NOT established.',
 'materials':{},'all_diagnostic_rows':len(rows),'input_hashes_equal':True}
phys=[]
for mat,nb in [('MoSSe',112),('InSe',160)]:
 rr=[r for r in rows if r['material']==mat and r['nb']==nb]
 selected=next(r for r in rr if r['matching']=='anchor' and r['b_eV']==.0000625)
 slopes=np.array([r['total'] for r in rr]);energy=np.array([r['energy'] for r in rr])
 cm=np.array(selected['Cminus']);cp=np.array(selected['Cplus']);mid=(cm+cp)/2
 relative=100*(abs(cp)-abs(cm))/abs(mid)
 summary['materials'][mat]={'nb':nb,'G_net_lower_upper':selected['total'],
 'G_range_min':slopes.min(0).tolist(),'G_range_max':slopes.max(0).tolist(),
 'energy_range_min':energy.min(0).tolist(),'energy_range_max':energy.max(0).tolist(),
 'Cminus':cm.tolist(),'Cplus':cp.tolist(),'relative_amplitude_change_percent_net_lower_upper':relative.tolist(),
 'upper_fraction_abs_layer_gate_changes':abs(selected['total'][2])/(abs(selected['total'][1])+abs(selected['total'][2])),
 'minimum_gap_eV':min(x[t]['gap_eV'] for r in rr for x in r['details'] for t in ('plus','minus')),
 'max_split_sum_error':float(max(np.max(abs(np.array(r['energy'])+r['matrix']-r['total'])) for r in rr)),
 'max_layer_sum_error':float(max(abs(r['total'][0]-sum(r['total'][1:])) for r in rr))}
 for i,layer in enumerate(('lower','upper'),1):
  phys.append(dict(material=mat,layer=layer,Cminus=cm[i],Cplus=cp[i],G=selected['total'][i],amplitude_change_percent=relative[i]))
with (ROOT/'layer_changes.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(phys[0]),lineterminator="\n");w.writeheader();w.writerows(phys)
summary['baseline_verification']=json.loads((ROOT/'baseline_verification.json').read_text())
assert summary['baseline_verification']['max_response_difference']<1e-5
assert max(m['max_split_sum_error'] for m in summary['materials'].values())<1e-6
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
