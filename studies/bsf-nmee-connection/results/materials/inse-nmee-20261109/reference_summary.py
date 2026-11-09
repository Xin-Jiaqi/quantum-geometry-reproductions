"""Keep both electrostatic taper definitions and all KS cutoffs in comparisons."""
from pathlib import Path
import json, csv
p = Path(__file__).resolve().parent
d = json.loads((p/'reference.json').read_text()); assert d['complete']
table = []
keys = sorted({(r['nb'], round(r['taper_end_A'], 6)) for r in d['points']['A0']['rows']})
for nb, end in keys:
    vals = {}
    for case, data in d['points'].items():
        vals[case] = {r['layer']: r['C_layer'] for r in data['rows'] if r['nb'] == nb and round(r['taper_end_A'], 6) == end}
        assert set(vals[case]) == {'lower', 'upper'}
    gl = (vals['Ap']['lower']-vals['Am']['lower'])/.2
    gu = (vals['Ap']['upper']-vals['Am']['upper'])/.2
    table.append(dict(nb=nb, taper_end_A=end, C0_lower=vals['A0']['lower'], C0_upper=vals['A0']['upper'],
                      C0_net=sum(vals['A0'].values()), G_lower=gl, G_upper=gu, G_net=gl+gu))
with (p/'reference_slopes.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(table[0]), lineterminator='\n'); w.writeheader(); w.writerows(table)
summary = {'rows': table, 'gate_potential_checks': d['gate_potential_checks'],
           'minimum_reference_gap_eV': min(r['gap_eV'] for case in d['points'].values() for r in case['rows']),
           'scope': d['scope']}
(p/'reference_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
