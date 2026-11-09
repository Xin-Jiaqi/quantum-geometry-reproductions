"""Summarize completed small channel records locally; no remote wave processing."""
from pathlib import Path
import json, csv, itertools
import numpy as np
p = Path(__file__).resolve().parent
d = json.loads((p/'channels.json').read_text()); assert d['complete']
rows = d['rows']; keys = sorted({(r['nb'], r['frozen_empty'], r['channel'], r['alpha']) for r in rows})
table = []
for nb, frozen, channel, alpha in keys:
    states = {r['case']: r for r in rows if (r['nb'], r['frozen_empty'], r['channel'], r['alpha']) == (nb, frozen, channel, alpha)}
    g = (np.array(states['Ap']['C_total_lower_upper'])-states['Am']['C_total_lower_upper'])/.2
    c = states['zero']['C_total_lower_upper']
    table.append(dict(nb=nb, frozen_empty=frozen, channel=channel, alpha=alpha,
                      C0_net=c[0], C0_lower=c[1], C0_upper=c[2], G_net=g[0], G_lower=g[1], G_upper=g[2],
                      upper_absolute_share=abs(g[2])/(abs(g[1])+abs(g[2])) if abs(g[1])+abs(g[2]) else None))
with (p/'channel_slopes.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(table[0]), lineterminator='\n'); w.writeheader(); w.writerows(table)

k = np.array(d['k']); cell = np.array(d['cell']); reciprocal = 2*np.pi*np.linalg.inv(cell).T
dist = np.min([np.linalg.norm((k+[x, y, 0])@reciprocal, axis=1) for x, y in itertools.product((-1, 0, 1), repeat=2)], axis=0)
gamma_k = np.linalg.norm(np.array([1/3, 1/3, 0])@reciprocal)
central = [r for r in rows if r['nb'] == 160 and r['frozen_empty'] == 8 and 'C_k_total_lower_upper' in r]
def maps(channel, alpha):
    rr = {r['case']: np.array(r['C_k_total_lower_upper']) for r in central if r['channel'] == channel and r['alpha'] == alpha}
    return {'C0': rr['zero'], 'G': (rr['Ap']-rr['Am'])/.2}
full = maps('full', 1.); novc = maps('VC', 0.)
sets = {'full': full, 'no_selected_VC': novc,
        'change_when_VC_restored': {name: full[name]-novc[name] for name in full}}
regions = []
for label, quantities in sets.items():
    for quantity, values in quantities.items():
        for radius in (.25, .40, .60):
            mask = dist < radius*gamma_k
            regions.append({'control': label, 'quantity': quantity, 'Gamma_radius_over_Gamma_K': radius,
                            'n_k': int(mask.sum()), 'total_k': len(k),
                            'region_signed_integral': (values[mask].sum(axis=0)/len(k)).tolist(),
                            'whole_signed_integral': values.mean(axis=0).tolist(),
                            'region_absolute_integrand_share': (abs(values[mask]).sum(axis=0)/abs(values).sum(axis=0)).tolist()})
out = {'order': 'total, lower, upper', 'rows': table, 'regions': regions,
       'scope': 'Uniform k-weight sum. Absolute-integrand share differs from signed response; channel intervention is not an additive contribution.',
       'minimum_gap_eV': min(r['gap_eV'] for r in rows),
       'max_half_probe_change': max(r.get('half_probe_max_change', 0.) for r in rows),
       'max_transport_singular_deviation': max(r['representation_max_singular_deviation'] for r in d['checks'])}
(p/'channel_summary.json').write_text(json.dumps(out, indent=2))
print('Central-window channels:')
for r in table:
    if r['nb'] == 160 and r['frozen_empty'] == 8:
        print(r)
print('Regional summaries saved; inspect signs and absolute shares separately.')
