"""Fixed zero-field valence/empty interlayer controls in the existing KM model."""
import argparse
import json
from pathlib import Path

import numpy as np

from compare_reversal_mechanism import (AuditBilayer, AuditSpec, ab_ba_shift,
                                       integration_grid, km_response)
from bsf_nmee.mixing_channels import split_channels


def run(grid, step):
    model = AuditBilayer(AuditSpec(shift=ab_ba_shift(.25), staggered_eV=.045,
                                  tunneling_eV=.020))
    direction = np.array([.5, np.sqrt(3) / 2])
    spins = [sum(direction[j] * model.spins[name + '_' + ax]
                 for j, ax in enumerate(['x', 'y']))
             for name in ['bottom', 'top', 'total']]
    k, weights, _ = integration_grid(model, *grid)
    subsets = ['none', 'VV', 'CC', 'VC', 'VV_CC', 'VV_VC', 'CC_VC', 'all']
    rows = {}
    for t in [10, 20]:
        for name in subsets:
            rows[f'{t}_{name}'] = {'C': np.zeros(3), 'vmax': -np.inf, 'cmin': np.inf}
    errors = []
    for start in range(0, len(k), 512):
        ks, w = k[start:start + 512], weights[start:start + 512]
        h = model.matrices(ks)
        v = sum(direction[j] * model.matrices(ks, axis=j) for j in range(2))
        intra, channels, _, check, vi, dv = split_channels(
            h, [np.arange(4), np.arange(4, 8)], 2, velocity=v)
        errors.append(check['reconstruction_eV'])
        for t in [10, 20]:
            for name in subsets:
                use = list(channels) if name == 'all' else ([] if name == 'none' else name.split('_'))
                hc = intra + t / 20 * sum((channels[c] for c in use), np.zeros_like(h))
                vc = vi + t / 20 * sum((dv[c] for c in use), np.zeros_like(v))
                val, _, vmax, cmin, _, _ = km_response(
                    hc, vc, spins, model.spins['total_z'], w, step, model.mono.g_factor)
                r = rows[f'{t}_{name}']
                r['C'] += val
                r['vmax'], r['cmin'] = max(r['vmax'], vmax), min(r['cmin'], cmin)
    for r in rows.values():
        r['gap_eV'] = r.pop('cmin') - r.pop('vmax')
        r['layer_sum_error'] = float(r['C'][0] + r['C'][1] - r['C'][2])
        r['C'] = r['C'].tolist()
    return {'system': 'KM', 'fraction': .25, 'mass_meV': 45, 'grid': grid,
            'nk': len(k), 'step_T': step, 'probe_angle_deg': 60,
            'units': '(nm V T)^-1', 'rows': rows,
            'reconstruction_eV': max(errors),
            'definition': 'Frozen zero-field isolated layer projectors; exact k derivatives of all channels. Hamiltonian interventions, not additive response contributions.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--grid', nargs=3, type=int, default=[41, 40, 64])
    parser.add_argument('--step', type=float, default=.002)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    data = run(args.grid, args.step)
    args.out.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps({name: r['C'][0] for name, r in data['rows'].items()}, indent=2))
