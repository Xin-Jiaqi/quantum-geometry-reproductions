"""Independent EB solver at resolvable broadening; do not trust smallest eta."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift, half_filled_spec
from bsf_nmee.static_observables import observables
from nmee_lab.ptbackend.liouville_general import eb_single


def main():
    k, w, meta = integration_grid(AuditBilayer(), 41, 45, 72)
    out = {'grid': meta, 'scope': 'equal adiabatic rates; numerical cross-check, not a relaxation model', 'models': {}}
    for name, shift in [('AB', ab_ba_shift(0)), ('generic_small', (.023, .041))]:
        spec = half_filled_spec(AuditSpec(shift=shift, tunneling_eV=.005), k)
        model = AuditBilayer(spec)
        row = {'static': observables(model, k, w, step_T=.0625)['C'][0][0], 'samples': []}
        for eta in (1e-3, 3e-4, 1e-4, 3e-5):
            value = 0j
            for start in range(0, len(k), 512):
                ks = k[start:start+512]
                v = eb_single(model.matrices(ks), model.matrices(ks, axis=0),
                              model.spins['total_z'], model.spins['total_x'],
                              spec.mu_eV, 0., model.mono.g_factor, 1j*eta, 1j*eta)
                value += w[start:start+512]@v
            row['samples'].append({'eta_eV': eta, 'real': float(value.real), 'imag': float(value.imag)})
        # The two smallest stable samples estimate the O(eta^2) intercept.
        a, b = row['samples'][-2:]
        row['eta_squared_intercept'] = (a['real']*b['eta_eV']**2-b['real']*a['eta_eV']**2)/(b['eta_eV']**2-a['eta_eV']**2)
        out['models'][name] = row
        print(name, row, flush=True)
    (ROOT/'results/bilayer_stage9_limit_check.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
