"""Validate archived evidence and the integrity of this discussion package."""
from pathlib import Path
import csv
import json
import re
import numpy as np

R = Path(__file__).resolve().parent


def rows(name):
    with (R/name).open() as f:
        return [{k: (v if k == 'grid' else float(v)) for k, v in r.items()}
                for r in csv.DictReader(f)]


fine = rows('normal_main.csv') + rows('gate_extra_main.csv')
near = rows('neighbourhood_robust.csv') + rows('gate_extra_robust.csv')
for dataset in [fine, near]:
    keys = [(r['lambda'], r['mass_meV'], r['t_meV'], r['U_top_minus_bottom_meV']) for r in dataset]
    assert len(keys) == len(set(keys)), 'duplicate control point in plotted data'
    assert len({r['grid'] for r in dataset}) == 1, 'mixed integration settings in one curve'
    for r in dataset:
        assert r['common_gap_meV'] > 0
        assert abs(r['bottom'] + r['top'] - r['total']) < 1e-12


def point(dataset, t, u, lam=.25, mass=45):
    return next(r for r in dataset if r['t_meV'] == t and r['U_top_minus_bottom_meV'] == u
                and r['lambda'] == lam and r['mass_meV'] == mass)


assert point(fine, 0, 0)['bottom'] < 0 < point(fine, 20, 0)['bottom']
assert point(fine, 16, 0)['bottom'] < 0
gate = sorted([r for r in fine if r['t_meV'] == 10 and r['U_top_minus_bottom_meV'] <= 5],
              key=lambda r: r['U_top_minus_bottom_meV'])
assert all(r['bottom'] < 0 < r['top'] for r in gate)
a, b = point(fine, 10, 4.5), point(fine, 10, 4.75)
assert a['total'] > 0 > b['total']
root = 4.5 + .25*a['total']/(a['total']-b['total'])
neighbour_report = []
for lam, mass in [(.24,45), (.25,45), (.26,45), (.25,44.9), (.25,45.1)]:
    assert point(near, 10, 0, lam, mass)['bottom'] < 0 < point(near, 20, 0, lam, mass)['bottom']
    assert point(near, 10, 4, lam, mass)['total'] > 0 > point(near, 10, 5, lam, mass)['total']
    subset = [r for r in near if r['lambda']==lam and r['mass_meV']==mass and r['t_meV']==10]
    assert all(r['bottom'] < 0 < r['top'] for r in subset)
    neighbour_report.append({'lambda': lam, 'mass_meV': mass, 'both_routes': True})

other = rows('normal_robust.csv')
for lam, mass in [(.30,45), (.25,46)]:
    assert all(point(other,t,0,lam,mass)['bottom'] < 0 for t in [10,20,30])

audit = json.loads((R/'forbidden_at_compensation.json').read_text())
for key in ['transverse','out_of_plane']:
    assert max(abs(audit['response'][key][k]) for k in ['bottom','top','total']) < 1e-10
assert max(row['max_error_eV'] for row in audit['mirror']) < 1e-12
# Independent original tensor must obey the transformed mirror selection rule.
rot = np.array([[.5,np.sqrt(3)/2],[-np.sqrt(3)/2,.5]])
offdiag = []
independent_errors = []
for case in json.loads((R/'tensor_check.json').read_text()):
    for layer, tensor in case['tensor'].items():
        c = rot@np.array(tensor)@rot.T
        offdiag.extend([abs(c[0,1]), abs(c[1,0])])
        if case['lambda'] == .25 and case['U_meV'] in [0, 12]:
            target = point(other, 10, case['U_meV'])[layer]
            independent_errors.append(abs(c[0,0] - target))
assert max(offdiag) < 1e-10
assert max(independent_errors) < 1e-10

mds = ['正文讨论稿-修订版.md', '理论定义与公式核对.md', '补充讨论与图.md', '批注回应与机制问题.md']
for name in mds:
    content = (R/name).read_text()
    assert content.count('$$') % 2 == 0
    assert not any(ord(c)<32 and c not in '\n\r\t' for c in content)
    for link in re.findall(r'\]\(([^)]+)\)', content):
        if '://' not in link and not link.startswith('#'):
            assert (R/link).exists(), (name, link)
    for block in re.findall(r'(?:^\|.*\n)+', content, flags=re.M):
        widths = {len(line.split('|')) for line in block.splitlines()}
        assert len(widths) == 1, (name, widths)

for name in ['Fig1_KM','FigS4_KM_AB_BA','FigS6_probe_and_range']:
    for ext in ['png','svg','pdf']:
        assert (R/'figures-20261106c'/f'{name}.{ext}').is_file()
    assert '<text' in (R/'figures-20261106c'/f'{name}.svg').read_text()

report = {'status':'pass', 'fine_points':len(fine), 'neighbour_points':len(near),
          'first_zero_bracket_meV':[4.5,4.75], 'linear_interpolation_meV':root,
          'min_fine_common_gap_meV':min(r['common_gap_meV'] for r in fine),
          'max_fine_layer_sum_error':max(abs(r['layer_sum_error']) for r in fine),
          'max_rotated_tensor_forbidden':max(offdiag), 'neighbours':neighbour_report,
          'max_independent_scalar_tensor_error':max(independent_errors),
          'negative_controls_retained':True, 'relative_links_and_tables':'pass'}
(R/'revision_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
