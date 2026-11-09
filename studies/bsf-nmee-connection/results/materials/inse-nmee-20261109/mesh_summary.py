"""Compare collected zero-gate mesh records; this is not a gate-slope audit."""
from pathlib import Path
import itertools
import json
import numpy as np

p = Path(__file__).resolve().parent
six = json.loads((p/'A0/result.json').read_text())
eight = json.loads((p/'A0_mesh8/result.json').read_text())
assert six['complete'] and eight['complete']
rows = []
for nb in ('128', '144', '160'):
    a = np.array(six['response'][nb]['standard']['C_total_lower_upper'])
    b = np.array(eight['response'][nb]['standard']['C_total_lower_upper'])
    rows.append({'bands': int(nb), 'C_6x6': a.tolist(), 'C_8x8': b.tolist(),
                 'relative_change': ((b-a)/abs(a)).tolist(),
                 'same_signs': (np.sign(a) == np.sign(b)).tolist(),
                 'half_probe_max_change_8x8': float(np.max(abs(b-eight['response'][nb]['half_b']['C_total_lower_upper'])))})
k = np.array(eight['k']); cell = np.array(eight['cell'])
reciprocal = 2*np.pi*np.linalg.inv(cell).T
dist = np.min([np.linalg.norm((k+[x,y,0])@reciprocal, axis=1)
               for x,y in itertools.product((-1,0,1), repeat=2)], axis=0)
gamma_k = np.linalg.norm(np.array([1/3,1/3,0])@reciprocal)
values = np.array(eight['C_k_total_lower_upper'])
np.testing.assert_allclose(values.mean(axis=0), rows[-1]['C_8x8'], rtol=1e-8, atol=1e-7)
regions = []
for radius in (.25, .4, .6):
    mask = dist < radius*gamma_k
    regions.append({'radius_over_Gamma_K': radius, 'n_k': int(mask.sum()), 'total_k': len(k),
                    'signed_integral': (values[mask].sum(axis=0)/len(k)).tolist(),
                    'absolute_integrand_share': (abs(values[mask]).sum(axis=0)/abs(values).sum(axis=0)).tolist()})
out = {'scope': 'Independent zero-gate SCF at 6x6 and 8x8; no gate sensitivity or full mesh convergence claim.',
       'order': 'total, lower, upper', 'gap_6x6_eV': six['export']['gap_eV'],
       'gap_8x8_eV': eight['export']['gap_eV'], 'rows': rows, 'regions_8x8': regions}
(p/'mesh_summary.json').write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
