"""Regression of the post-relaxation guard on actual saved structures."""
from pathlib import Path
import json
from ase.io import read
from geometry_check import validate_fixed_registry

root = Path(__file__).resolve().parent
out = {}
for label in ('e000', 'em02', 'ep02'):
    initial = read(root / label / 'relax/POSCAR')
    relaxed = read(root / label / 'relax/CONTCAR')
    out[label] = validate_fixed_registry(initial, relaxed)
    out[label]['actual_all_ions_free'] = not bool(initial.constraints)

# A real 0.01-A registry change must not be accepted as numerical roundoff.
shifted = relaxed.copy()
shifted.positions[0, 0] += .01
try:
    validate_fixed_registry(initial, shifted)
except ValueError:
    out['reject_0p01_A_displacement'] = True
else:
    raise AssertionError('Guard accepted a physically shifted atom')

# Equivalent periodic coordinates represent the same stacking.
wrapped = relaxed.copy()
wrapped.positions[0] += wrapped.cell[0]
out['periodic_image_accepted'] = validate_fixed_registry(initial, wrapped)
(root / 'registry_validation.json').write_text(json.dumps(out, indent=2) + '\n')
print('Registry guard accepts roundoff/periodic images and rejects a real displacement')
