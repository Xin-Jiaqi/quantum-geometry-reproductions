"""Check completed inputs and geometry without imposing the predicted outcome."""
from pathlib import Path
import json
import numpy as np
from ase.io import read
import spglib

ROOT = Path(__file__).resolve().parent
jobs = json.loads((ROOT / 'jobs.json').read_text())
states = json.loads((ROOT / 'final_job_status.json').read_text())
assert len(jobs) == 7 and len(states) == 7
assert all(v['job_state'] == 'F' and v['Exit_status'] == '0' for v in states.values())
reference = json.loads((ROOT / 'reference_status.json').read_text())
assert reference['job_state'] == 'F' and reference['Exit_status'] == '0'
base = read(ROOT / 'd000_fm032/seed/POSCAR')
lower = base.positions[:, 2] < base.cell[2, 2] / 2
checks = []
for case, job in jobs.items():
    path = ROOT / case
    meta = json.loads((path / 'input.json').read_text())
    export = json.loads((path / 'export.json').read_text())
    atoms = read(path / 'seed/POSCAR')
    expected = base.positions.copy()
    expected[:, 2] += np.where(lower, -meta['delta_d_A']/2, meta['delta_d_A']/2)
    assert np.max(abs(atoms.positions - expected)) < 1e-9
    assert np.max(abs(atoms.cell.array - base.cell.array)) < 1e-9
    assert np.array_equal(atoms.numbers, base.numbers)
    assert meta['job'] == job and meta['node'].startswith('csuper')
    assert export['nk'] == 36 and export['bands'] == 112 and export['nocc'] == 48
    assert export['metric_max'] < 1e-4 and export['link_reverse_error'] < 1e-6
    assert export['wave_link_gauge_error'] < .08
    inc = dict(line.split('=', 1) for line in (path / 'seed/INCAR').read_text().splitlines() if '=' in line)
    assert inc['BEXT'].strip() == '0 0 0'
    assert np.isclose(float(inc['EFIELD']), -meta['physical_Fz_Vnm']/10)
    assert inc['NSW'].strip() == '0' and inc['ICHARG'].strip() == '1'
    sym = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-5)
    assert sym.pointgroup == '3m'
    checks.append(dict(case=case, geometric_point_group=sym.pointgroup,
        rigid_layer_compression_verified=True, zero_B_reference=True,
        input_job=job, export_checks_passed=True))
(ROOT / 'validation.json').write_text(json.dumps({'checks':checks,
    'scope':'Input and numerical integrity only; the physical prediction is not a pass/fail condition.'}, indent=2))
print('Seven new cases, geometry, control fields, compute allocations, exports and final job exits verified.')
