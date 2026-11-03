"""Vary hybridization at one fixed point of the physical AB-to-BA segment."""
import json
from dataclasses import asdict
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift, half_filled_spec
from bsf_nmee.static_observables import observables


def main():
    k, w, meta = integration_grid(AuditBilayer(), 61, 60, 96)
    result = {'grid': meta, 'lambda': .25,
              'scope': 'fixed stacking, hypothetical hopping control; no pressure calibration',
              'scan': [], 'fine': []}
    for t in (0., .0025, .005, .0075, .01, .0125, .015, .02, .025, .03):
        spec = half_filled_spec(AuditSpec(shift=ab_ba_shift(.25), tunneling_eV=t), k)
        result['scan'].append({'parameters': asdict(spec),
                              **observables(AuditBilayer(spec), k, w, step_T=.0625)})
        print('coupling', t, result['scan'][-1]['difference'], flush=True)
        (ROOT/'results/bilayer_stage9_coupling.json').write_text(json.dumps(result, indent=2)+'\n')
    k, w, meta_fine = integration_grid(AuditBilayer(), 101, 100, 144)
    result['fine_grid'] = meta_fine
    for t in (.005, .02):
        spec = half_filled_spec(AuditSpec(shift=ab_ba_shift(.25), tunneling_eV=t), k)
        result['fine'].append({'parameters': asdict(spec),
                              **observables(AuditBilayer(spec), k, w, step_T=.03125)})
    (ROOT/'results/bilayer_stage9_coupling.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
