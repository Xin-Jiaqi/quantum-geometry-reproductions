from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from ase.io import read, write
from night_checks import geometry_acceptance, prepare_inputs, response_acceptance, compare_meshes, digest

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent/'full-tensor-pilot-20261110'


def test_precision_artifact_does_not_relax_actual_xy_constraint(tmp_path):
    initial = PILOT/'inputs7268/POSCAR.A'
    final = PILOT/'pilot7268/A_zrelax/CONTCAR'
    record = json.loads((PILOT/'pilot7268/A_zrelax/result.json').read_text())
    assert not record['passed']
    assert geometry_acceptance(initial, final, record)['passed']
    moved = read(final)
    moved.positions[0, 0] += 2e-5
    candidate = tmp_path/'CONTCAR'
    write(candidate, moved, format='vasp', direct=True)
    record['CONTCAR_sha256'] = digest(candidate)
    checked = geometry_acceptance(initial, candidate, record)
    assert not checked['checks']['xy_fixed']
    assert not checked['passed']


def good_response():
    d = json.loads((PILOT/'pilot7268/A0/result.json').read_text())
    d.update(label='A', max_force_eV_A=.005, gauge_audit={'passed': True})
    assert response_acceptance(d, .03)['passed']
    return d


def test_closed_gap_and_nonfinite_half_probe_are_rejected():
    for gap in (0.0, float('nan')):
        d = good_response()
        d['responses']['112']['half_b']['gap_eV'] = gap
        assert not response_acceptance(d, .03)['passed']


def test_forbidden_spin_and_band_direction_disagreement_are_rejected():
    d = good_response()
    d['responses']['112']['standard']['C'][1][2][0] = 10.
    assert not response_acceptance(d)['passed']
    d = good_response()
    for probe in ('standard', 'half_b'):
        for i in (0, 1):
            d['responses']['112'][probe]['C'][1][i][i] *= -1
    assert not response_acceptance(d)['band_window_signs_agree']


def test_wrong_electron_count_is_rejected():
    d = good_response()
    d['export']['nocc'] = 50
    assert not response_acceptance(d)['passed']


def test_matched_isolates_preserve_internal_geometry_and_remain_distinct(tmp_path):
    original = read(PILOT/'pilot7268/A_zrelax/CONTCAR')
    meta = prepare_inputs(PILOT/'pilot7268/A_zrelax/CONTCAR', tmp_path/'inputs')
    assert meta['max_isolated_z_difference_A'] > .001
    for name, inds in [('lower', [0, 2, 4]), ('upper', [1, 3, 5])]:
        mono = read(tmp_path/'inputs'/f'POSCAR.{name}')
        delta = mono.positions-original.positions[inds]
        assert np.max(np.abs(delta-delta.mean(axis=0))) < 1e-12
        assert mono.get_chemical_symbols() == ['Mo', 'S', 'Se']


def test_equal_gate_slopes_do_not_force_a_sensitive_layer():
    observation = {}
    for label in ('A', 'lower', 'upper'):
        n = 3 if label == 'A' else 1
        c = np.zeros((n, 3, 2)); c[:, 0, 0] = 1; c[:, 1, 1] = 1
        row = {'C0': c.tolist(), 'G': c.tolist(), 'G_half': c.tolist()}
        observation[label] = {**row, 'by_band': {'64': deepcopy(row)}}
    result = compare_meshes(observation, deepcopy(observation))
    assert result['agreement']
    assert result['gate_sensitive_layer'] == 'unresolved'


def test_worker_refuses_execution_without_scheduler(tmp_path):
    env = dict(os.environ)
    env.pop('PBS_JOBID', None); env.pop('PBS_NODEFILE', None)
    proc = subprocess.run([sys.executable, str(HERE/'night.py')], cwd=tmp_path,
                          env=env, capture_output=True, text=True, timeout=10)
    assert proc.returncode != 0
    assert 'PBS compute allocation required' in proc.stderr
    assert not list(tmp_path.iterdir())


def test_band_window_gate_reversal_is_unresolved_even_when_meshes_identical():
    observation = {}
    for label in ('A', 'lower', 'upper'):
        n = 3 if label == 'A' else 1
        c = np.zeros((n, 3, 2)); c[:, 0, 0] = np.arange(1,n+1); c[:, 1, 1] = c[:, 0, 0]
        row = {'C0': c.tolist(), 'G': c.tolist(), 'G_half': c.tolist()}
        observation[label] = {**row, 'by_band': {'48': deepcopy(row), '64': deepcopy(row)}}
    observation['A']['by_band']['48']['G'][1][0][0] *= -1
    result = compare_meshes(observation, deepcopy(observation))
    assert result['agreement']
    assert not result['band_windows_adequate']
    assert result['gate_sensitive_layer'] == 'unresolved'
