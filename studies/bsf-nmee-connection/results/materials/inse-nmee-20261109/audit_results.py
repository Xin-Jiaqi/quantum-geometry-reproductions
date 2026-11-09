"""Audit the small archived results without opening remote wavefunctions."""
from pathlib import Path
import ast
import json
import numpy as np

p = Path(__file__).resolve().parent
names = ('mono0', 'monom', 'monop', 'A0', 'B0', 'Am', 'Ap')
data = {n: json.loads((p/n/'result.json').read_text()) for n in names}
assert all(d['complete'] and d['fixed_ions'] for d in data.values())
assert all(d['node'].startswith('csuper') for d in data.values())
assert [data[n]['F_external_V_nm'] for n in names] == [0, -.1, .1, 0, 0, -.1, .1]
sum_error = 0.; gaps = []
for d in data.values():
    for pair in d['response'].values():
        for rec in pair.values():
            c = rec['C_total_lower_upper']
            sum_error = max(sum_error, abs(c[0]-c[1]-c[2]))
            gaps.append(rec['gap_eV'])
assert sum_error < 1e-8 and min(gaps) > 0
channels = json.loads((p/'channels.json').read_text()); assert channels['complete']
restoration_error = 0.
for row in channels['rows']:
    if row['channel'] != 'full': continue
    name = 'A0' if row['case'] == 'zero' else row['case']
    expected = data[name]['response'][str(row['nb'])]['standard']['C_total_lower_upper']
    restoration_error = max(restoration_error, float(max(abs(np.array(row['C_total_lower_upper'])-expected))))
assert restoration_error < 1e-8
for script in p.glob('*.py'): ast.parse(script.read_text())
out = {'seven_points_complete': True, 'all_science_nodes_are_compute_nodes': True,
       'net_equals_layer_sum_max_error': sum_error, 'full_channel_restoration_max_error': restoration_error,
       'minimum_probe_gap_eV': min(gaps), 'python_sources_parse': True,
       'scope': 'Data/provenance consistency; not a proof of physical or mesh convergence.'}
(p/'audit_results.json').write_text(json.dumps(out, indent=2)); print(json.dumps(out, indent=2))
