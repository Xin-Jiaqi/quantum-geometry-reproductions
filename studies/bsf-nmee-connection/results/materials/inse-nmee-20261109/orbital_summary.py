"""Summarize collected Gamma character weights across representations."""
from pathlib import Path
import json

p = Path(__file__).resolve().parent
d = json.loads((p/'orbital_snapshot.json').read_text())
rows = []
for pair in ([69,70], [71,72], [73,74], [75,76], [77,78], [79,80]):
    subset = [r for r in d['rows'] if r['bands_1_based'] == pair]
    assert len(subset) == 9
    rows.append({'bands_1_based': pair, 'energy_eV': subset[0]['energy_eV'],
                 'character_ranges': {key: [min(r['character'][key] for r in subset),
                                             max(r['character'][key] for r in subset)]
                                      for key in subset[0]['character']},
                 'lower_weight_range': [min(r['lower_weight'] for r in subset),
                                        max(r['lower_weight'] for r in subset)]})
out = {'scope': d['scope'],
       'basis': 'Orthonormalized selected atomic trials; descriptive, not causal response decomposition.',
       'rows': rows}
(p/'orbital_summary.json').write_text(json.dumps(out, indent=2))
