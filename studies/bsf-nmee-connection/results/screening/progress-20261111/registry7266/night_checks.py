"""Explicit scientific gates for the single, bounded 7268 night allocation."""
import hashlib
import json
from pathlib import Path
import numpy as np
import spglib
from ase.io import read, write


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def geometry_acceptance(initial_path, contcar_path, relaxation_record):
    original, final = read(initial_path), read(contcar_path)
    delta = final.positions-original.positions
    delta[:, 2] -= delta[:, 2].mean()
    ds = spglib.get_symmetry_dataset((final.cell.array, final.get_scaled_positions(), final.numbers), symprec=1e-5)
    checks = {
        'source_hash': digest(contcar_path) == relaxation_record['CONTCAR_sha256'],
        'species_and_order': original.get_chemical_symbols() == final.get_chemical_symbols(),
        'vasp_completed': relaxation_record['returncode'] == 0 and relaxation_record['complete'],
        'ionic_convergence': relaxation_record['ionic_converged'],
        'all_electronic_steps_converged': relaxation_record['electronic_converged_steps'] >= relaxation_record['ionic_steps'],
        'ionic_budget': relaxation_record['ionic_steps'] <= 12,
        'force': relaxation_record['max_force_eV_A'] < .03,
        'xy_fixed': bool(np.max(np.abs(delta[:, :2])) < 1e-6),
        'cell_fixed': bool(np.max(np.abs(final.cell.array-original.cell.array)) < 1e-8),
        'z_bounded': bool(np.max(np.abs(delta[:, 2])) <= .15),
        'point_group': ds.pointgroup == '3m',
        'no_internal_layer_exchange': not any(r[2, 2] < 0 for r in ds.rotations),
    }
    return {'passed': all(checks.values()), 'checks': checks,
            'coordinate_source': 'high-precision CONTCAR; forces/convergence remain from OUTCAR',
            'original_wrapper_passed': relaxation_record['passed'],
            'coarse_OUTCAR_xy_drift_A': relaxation_record['max_xy_displacement_A'],
            'precise_xy_drift_A': float(np.max(np.abs(delta[:, :2]))),
            'precise_centered_z_motion_A': float(np.max(np.abs(delta[:, 2]))),
            'max_force_eV_A': relaxation_record['max_force_eV_A'],
            'point_group': ds.pointgroup, 'layer_exchanging_operations': 0,
            'source_sha256': {'initial': digest(initial_path), 'CONTCAR': digest(contcar_path)},
            'adjudication': 'OUTCAR printed-coordinate false rejection; original failed record retained; no tolerance relaxed'}


def prepare_inputs(contcar, destination):
    destination.mkdir(exist_ok=False)
    atoms = read(contcar)
    atoms.set_constraint()
    order = np.argsort(atoms.positions[:, 2])
    lower, upper = np.array_split(order, 2)
    midpoint = (atoms.positions[lower, 2].max()+atoms.positions[upper, 2].min())/2
    shift = atoms.cell[2, 2]/2-midpoint
    atoms.positions[:, 2] += shift
    write(destination/'POSCAR.A', atoms, format='vasp', direct=True, sort=False)
    monos = []
    for name, indices in [('lower', lower), ('upper', upper)]:
        mono = atoms[np.sort(indices)]
        mono.positions[:, 2] += mono.cell[2, 2]/2-mono.positions[:, 2].mean()
        write(destination/f'POSCAR.{name}', mono, format='vasp', direct=True, sort=False)
        monos.append(mono)
    difference = monos[1].positions[:, 2]-monos[0].positions[:, 2]
    return {'bilayer_common_z_shift_A': float(shift), 'layer_boundary_fractional_z': .5,
            'both_isolated_layers_computed': True,
            'isolated_centered_z_difference_A': difference.tolist(),
            'max_isolated_z_difference_A': float(np.max(np.abs(difference))),
            'preserved': 'cell, xy registry, species order, all internal atomic coordinates; only common translations for centering',
            'source_CONTCAR_sha256': digest(contcar),
            'input_hashes': {f.name: digest(f) for f in destination.iterdir()}}


def response_acceptance(record, force_limit=.1):
    bilayer = record['label'] == 'A'
    scopes = [0, 1, 2] if bilayer else [0]
    rows, signs = [], []
    for nb, values in record['responses'].items():
        standard, half = values['standard'], values['half_b']
        c = np.asarray(standard['C'])[scopes]
        h = np.asarray(half['C'])[scopes]
        scale = max(float(np.max(np.abs(c[:, 0, 0]))), 1e-8)
        finite = bool(np.isfinite(c).all() and np.isfinite(h).all()
                      and np.isfinite([standard['gap_eV'], half['gap_eV']]).all())
        forbidden = float(np.max(np.abs(c[:, [0, 1, 2, 2], [1, 0, 0, 1]])))
        isotropy = float(np.max(np.abs(c[:, 0, 0]-c[:, 1, 1])))
        probe = float(np.max(np.abs(c-h)))
        closure = float(standard['layer_sum_max_abs'])
        passed = (finite and min(standard['gap_eV'], half['gap_eV']) > .05
                  and forbidden < .02*scale and isotropy < .05*scale
                  and probe < .01*scale and closure < 1e-6*scale
                  and standard['operator_layer_sum_max_abs'] < 1e-7)
        rows.append({'nb': int(nb), 'passed': passed, 'Cxx': c[:, 0, 0].tolist(),
                     'scale': scale, 'forbidden': forbidden, 'isotropy': isotropy,
                     'probe_difference': probe, 'gap_eV': standard['gap_eV'], 'layer_closure': closure})
        signs.append(np.sign(c[:, 0, 0]))
    same_signs = bool(signs) and all(np.array_equal(signs[0], s) for s in signs[1:])
    force_pass = bool(np.isfinite(record['max_force_eV_A']) and record['max_force_eV_A'] <= force_limit)
    gauge_pass = bool(record['gauge_audit']['passed'])
    electron_count = record['export']['nocc'] == (48 if bilayer else 24)
    return {'passed': all(r['passed'] for r in rows) and same_signs and force_pass and gauge_pass and electron_count,
            'rows': rows, 'band_window_signs_agree': same_signs,
            'force_pass': force_pass, 'force_limit_eV_A': force_limit, 'gauge_pass': gauge_pass,
            'electron_count_matches_protocol': electron_count,
            'limits': 'Pilot numerical gates; not final absolute-amplitude convergence.'}


def mesh_observables(records):
    result = {}
    for label in ('A', 'lower', 'upper'):
        scope = slice(0, 3) if label == 'A' else slice(0, 1)
        bands = {}
        for nb in records[label+'0']['responses']:
            arrays = {s: np.array(records[label+s]['responses'][nb]['standard']['C'])[scope]
                      for s in ('0', 'm', 'p')}
            half = {s: np.array(records[label+s]['responses'][nb]['half_b']['C'])[scope]
                    for s in ('m', 'p')}
            bands[nb] = {'C0': arrays['0'].tolist(), 'G': ((arrays['p']-arrays['m'])/.2).tolist(),
                         'G_half': ((half['p']-half['m'])/.2).tolist()}
        highest = max(bands, key=int)
        result[label] = {**bands[highest], 'by_band': bands,
                         'G_definition': '[C(+0.10)-C(-0.10)]/0.20 V/nm'}
    return result


def band_diagnostics(observation):
    checks, rank_errors = [], []
    for label in ('A', 'lower', 'upper'):
        use = slice(1, 3) if label == 'A' else slice(0, 1)
        reference = observation[label]
        for nb, row in reference['by_band'].items():
            for quantity, tolerance in [('C0', .05), ('G', .10)]:
                a, b = np.array(row[quantity])[use], np.array(reference[quantity])[use]
                scale = max(float(np.max(np.abs(b[:, 0, 0]))), 1e-8)
                difference = float(np.max(np.abs(a-b)))
                checks.append(difference <= tolerance*scale)
            g, gh = np.array(row['G'])[use, 0, 0], np.array(row['G_half'])[use, 0, 0]
            refg = np.array(reference['G'])[use, 0, 0]
            floor = max(float(np.max(np.abs(g-gh))), 1e-8)
            # Slopes below probe uncertainty are unresolved, not assigned a sign.
            checks.append(bool(np.all(np.abs(g) > 2*floor) and np.array_equal(np.sign(g), np.sign(refg))))
            if label == 'A':
                rank_errors.append(max(float(np.max(np.abs(g-refg))), floor))
                checks.append(bool(np.argmax(np.abs(g)) == np.argmax(np.abs(refg))))
    return {'adequate': all(checks), 'observed_band_or_probe_change': max(rank_errors, default=1e-8)}


def compare_meshes(previous, current):
    checks = []
    for label in ('A', 'lower', 'upper'):
        for quantity, tolerance in [('C0', .05), ('G', .10)]:
            old, new = np.array(previous[label][quantity]), np.array(current[label][quantity])
            use = slice(1, 3) if label == 'A' else slice(0, 1)
            scale = max(float(np.max(np.abs(new[use, 0, 0]))), 1e-8)
            difference = float(np.max(np.abs(new[use]-old[use])))
            checks.append({'label': label, 'quantity': quantity, 'absolute_difference': difference,
                           'scale': scale, 'relative_to_largest_physical_layer': difference/scale,
                           'tolerance': tolerance, 'agreement': difference <= tolerance*scale})
    old_band, new_band = band_diagnostics(previous), band_diagnostics(current)
    band_ok = old_band['adequate'] and new_band['adequate']
    g = np.array(current['A']['G'])[1:3, 0, 0]
    old_g = np.array(previous['A']['G'])[1:3, 0, 0]
    uncertainty = float(np.max(np.abs(g-old_g)))
    band_uncertainty = max(old_band['observed_band_or_probe_change'], new_band['observed_band_or_probe_change'])
    margin = float(abs(abs(g[0])-abs(g[1])))
    resolved = bool(band_ok and margin > 2*max(uncertainty, band_uncertainty)
                    and np.argmax(np.abs(g)) == np.argmax(np.abs(old_g)))
    return {'agreement': all(c['agreement'] for c in checks), 'checks': checks,
            'band_windows_adequate': band_ok, 'band_diagnostics': {'previous': old_band, 'current': new_band},
            'gate_sensitive_layer': ('lower' if abs(g[0]) > abs(g[1]) else 'upper') if resolved else 'unresolved',
            'rank_resolved_against_observed_mesh_change': resolved,
            'rank_margin': margin, 'observed_mesh_change': uncertainty,
            'limits': 'Observed mesh, band and probe differences are diagnostics, not statistical errors; net residual needs separate review.'}
