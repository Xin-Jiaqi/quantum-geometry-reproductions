"""Audit collected small inputs/results locally; no remote scientific processing."""
from pathlib import Path
import json, hashlib
import numpy as np
import spglib
from ase.io import read
from geometry_check import validate_fixed_registry

root = Path(__file__).resolve().parent
records = {}
zero = read(root / 'e000/relax/POSCAR')
for label, strain in [('em02', -.02), ('e000', 0.), ('ep02', .02)]:
    base = root / label
    result = json.loads((base / 'result.json').read_text())
    assert result['stage'] == 'complete'
    initial = read(base / 'relax/POSCAR')
    relaxed = read(base / 'relax/CONTCAR')
    assert not initial.constraints, 'Actual inputs allowed all ionic coordinates'
    expected = zero.cell.array.copy()
    expected[:2] *= 1 + strain
    assert np.max(abs(initial.cell.array - expected)) < 1e-7
    registry = validate_fixed_registry(initial, relaxed)
    symmetry = spglib.get_symmetry_dataset(
        (relaxed.cell.array, relaxed.get_scaled_positions(), relaxed.numbers), symprec=1e-5)
    assert symmetry.pointgroup == '3m'
    bilayer = read(base / 'A/seed/POSCAR')
    isolated = read(base / 'mono/seed/POSCAR')
    lower = bilayer[bilayer.positions[:, 2] < bilayer.cell[2, 2]/2]
    lower.positions[:, 2] -= lower.positions[0, 2]
    isolated.positions[:, 2] -= isolated.positions[0, 2]
    assert lower.get_chemical_symbols() == isolated.get_chemical_symbols()
    geometry_error = float(np.max(abs(lower.positions - isolated.positions)))
    assert geometry_error < 1e-6
    rec = {'strain': strain, 'pointgroup': symmetry.pointgroup,
           'spacegroup_supercell': symmetry.international, 'registry': registry,
           'matched_isolated_geometry_error_A': geometry_error,
           'max_force_components_eVA': {k: v for k, v in result['geometry'].items() if 'force' in k},
           'responses': {}}
    for material in ('A', 'mono'):
        seed = base / material / 'seed'
        finished = json.loads((seed / 'finished.json').read_text())
        for filename, expected_hash in finished['input_sha256'].items():
            assert hashlib.sha256((seed / filename).read_bytes()).hexdigest() == expected_hash
        incar = (base / material / 'seed/INCAR').read_text()
        assert 'BEXT=0 0 0' in incar and 'EFIELD=0\n' in incar
        ex = json.loads((base / material / 'export.json').read_text())
        assert ex['gap_eV'] > .05 and ex['link_reverse_error'] < 1e-6
        assert ex['metric_max'] < 1e-4 and ex['wave_link_gauge_error'] < .08
        rec[material + '_zeroB_gap_eV'] = ex['gap_eV']
        rec[material + '_input_hashes_match_remote'] = True
        rr = result[material + '_response']
        rec['responses'][material] = {}
        for nb, data in rr.items():
            c = np.array(data['standard']['C_total_lower_upper'])
            half = np.array(data['half_b']['C_total_lower_upper'])
            assert abs(c[0] - c[1] - c[2]) < 1e-8
            relevant = slice(None) if material == 'A' else slice(0, 1)
            rec['responses'][material][nb] = {
                'C': c.tolist() if material == 'A' else [float(c[0])],
                'max_relevant_half_step_change': float(abs(c[relevant] - half[relevant]).max())}
    records[label] = rec
(root / 'verification.json').write_text(json.dumps(records, indent=2) + '\n')
print('Verified same stacking, matched layer geometry, insulation and response consistency in all three cases')
