"""Collect actual response points; bracket and interpolate zeros without extrapolation."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def zeros(rows):
    """Piecewise-linear estimates are descriptive, not extra calculated points."""
    rows = sorted(rows, key=lambda r: r['Fz_Vnm'])
    answer = []
    for a, b in zip(rows[:-1], rows[1:]):
        ya, yb = a['C_total'], b['C_total']
        if ya * yb < 0:
            x = a['Fz_Vnm'] - ya * (b['Fz_Vnm'] - a['Fz_Vnm']) / (yb - ya)
            answer.append({'bracket_Vnm': [a['Fz_Vnm'], b['Fz_Vnm']],
                           'interpolated_Vnm': x})
    return answer


def main():
    cases = json.loads((ROOT / 'reference_cases.json').read_text())
    full, differences, export_checks = [], [], []
    for case in cases:
        name = case['case']
        path = ROOT / ('baseline.json' if name == 'd000_f000' else name + '/result.json')
        result = json.loads(path.read_text())
        for cutoff, response in result['results'].items():
            standard = response['standard']
            total, lower, upper = standard['C_total_lower_upper']
            assert np.isclose(total, lower + upper, atol=1e-9)
            assert standard['gap_eV'] > 0
            full.append(dict(case=name, delta_d_A=case['delta_d_A'], Fz_Vnm=case['Fz_Vnm'],
                cutoff=int(cutoff), C_total=total, C_lower=lower, C_upper=upper,
                gap_min_eV=standard['gap_eV']))
            if 'half_b' in response:
                differences.extend(abs(np.array(standard['C_total_lower_upper']) -
                                       response['half_b']['C_total_lower_upper']))
        inputs = ROOT / name / 'input.json'
        if inputs.exists():
            metadata = json.loads(inputs.read_text())
            for filename, digest in metadata['input_sha256'].items():
                assert hashlib.sha256((ROOT / name / 'seed' / filename).read_bytes()).hexdigest() == digest
            export_checks.append(json.loads((ROOT / name / 'export.json').read_text()))
    static = json.loads((ROOT / 'electrostatic_reference.json').read_text())['rows']
    assert len(full) == 18 and len(static) == 72
    for rows, filename in [(full, 'full_response.csv'), (static, 'electrostatic_response.csv')]:
        with (ROOT / filename).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    root_records = []
    for label, rows, options in [('full', full, [(96, None), (112, None)]),
                                 ('electrostatic', static, [(n, t) for n in (40, 48, 56, 64) for t in (3.2, 3.8)])]:
        for cutoff, taper in options:
            for distance in (0, -.1, -.2):
                selected = [r for r in rows if r['cutoff'] == cutoff and r['delta_d_A'] == distance
                            and (taper is None or r['taper_end_A'] == taper)]
                assert len(selected) == 3
                endpoints = [r for r in selected if r['Fz_Vnm'] in (0, -.5)]
                root_records.append(dict(model=label, cutoff=cutoff, taper_end_A=taper, delta_d_A=distance,
                    three_point_bracket=zeros(selected), endpoints_only=zeros(endpoints),
                    layer_signs_at_actual_points=[dict(Fz_Vnm=r['Fz_Vnm'], lower=int(np.sign(r['C_lower'])),
                         upper=int(np.sign(r['C_upper']))) for r in sorted(selected, key=lambda r:r['Fz_Vnm'])]))
    summary = {'units': '1e-9 (nm V T)^-1, Pauli spin density per area',
        'scope': 'Fixed-registry frozen-ion, static spin-only full-KS probe response; gate self-consistent.',
        'root_scope': 'Three actual fields per geometry. Interpolations are estimates; cutoff/taper spread is not statistical uncertainty.',
        'full_rows': full, 'compensation': root_records,
        'checks': {'min_full_gap_eV': min(r['gap_min_eV'] for r in full),
                   'min_reference_gap_eV': min(r['gap_min_eV'] for r in static),
                   'max_half_b_change': float(max(differences)),
                   'all_full_lower_negative': all(r['C_lower'] < 0 for r in full),
                   'all_full_upper_positive': all(r['C_upper'] > 0 for r in full),
                   'new_input_hash_sets_verified': len(export_checks)}}
    summary['zero_gate_and_gate_leverage'] = []
    for cutoff in (96, 112):
        for distance in (0, -.1, -.2):
            selected = [r for r in full if r['cutoff'] == cutoff and r['delta_d_A'] == distance]
            at_zero = next(r for r in selected if r['Fz_Vnm'] == 0)
            at_negative = next(r for r in selected if r['Fz_Vnm'] == -.5)
            summary['zero_gate_and_gate_leverage'].append(dict(cutoff=cutoff, delta_d_A=distance,
                zero_gate_total=at_zero['C_total'],
                gate_secant_slope=(at_negative['C_total'] - at_zero['C_total']) / -.5,
                lower_change=at_negative['C_lower'] - at_zero['C_lower'],
                upper_change=at_negative['C_upper'] - at_zero['C_upper']))
    (ROOT / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary['checks'], indent=2))
    for model in ('full', 'electrostatic'):
        for distance in (0, -.1, -.2):
            vals = [z['interpolated_Vnm'] for row in root_records if row['model'] == model and row['delta_d_A'] == distance
                    for z in row['three_point_bracket']]
            print(model, distance, 'root range', (min(vals), max(vals)) if vals else 'no bracket')


if __name__ == '__main__':
    main()
