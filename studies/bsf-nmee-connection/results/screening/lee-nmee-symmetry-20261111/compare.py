"""Project the in-plane current-driven and E times Bz spin tensors for actual BiDB cells."""
from pathlib import Path
import json
import numpy as np
import spglib
from ase.io import read

here = Path(__file__).resolve().parent
structures = here.parent / "bidb-pilot-plan-20261110" / "structures"
transverse = np.array([[0., -1.], [1., 0.]])
longitudinal = np.eye(2)


def project(path, layers, magnetic):
    atoms = read(path)
    lattice = atoms.cell.array.T
    operations = spglib.get_symmetry_dataset(
        (atoms.cell, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-5)
    size = 4 * layers
    projector = np.zeros((size, size))
    for fractional in operations.rotations:
        rotation = lattice @ fractional @ np.linalg.inv(lattice)
        inplane = rotation[:2, :2]
        spin = np.linalg.det(rotation) * inplane
        bz = np.linalg.det(rotation) * rotation[2, 2]
        for i, basis in enumerate(np.eye(size)):
            tensor = basis.reshape(layers, 2, 2)
            image = spin @ tensor @ inplane.T
            if layers == 2 and rotation[2, 2] < 0:
                image = image[::-1]
            projector[:, i] += (image * (bz if magnetic else 1.)).ravel()
    projector /= len(operations.rotations)
    polarization = longitudinal if magnetic else transverse
    # Verify the whole projected space, including the interlayer relation.
    trial = []
    if operations.pointgroup == '3m':
        for layer in range(layers):
            tensor = np.zeros((layers, 2, 2))
            tensor[layer] = polarization
            trial.append(tensor.ravel())
    elif operations.pointgroup == '-3m' and layers == 2:
        trial.append(np.stack([polarization, -polarization]).ravel())
    expected = sum((np.outer(v, v) / np.dot(v, v) for v in trial), np.zeros_like(projector))
    assert np.allclose(expected, projector, atol=1e-7)
    return {'point_group': operations.pointgroup, 'allowed_coefficients': len(trial),
            'in_plane_spin_for_Ex': 'forbidden' if not trial else ('along x' if magnetic else 'along y'),
            'layers': ('independent by symmetry' if len(trial) == 2 else
                       'equal and opposite' if layers == 2 and trial else 'single layer')}


results = {}
for name in ('mono-1929', 'mono-2652', 'mono-7265', 'bilayer-1930', 'bilayer-2654',
             'bilayer-1931', 'bilayer-7269', 'bilayer-7268', 'bilayer-864', 'bilayer-741'):
    path = structures / (name + '.vasp')
    layers = 2 if name.startswith('bilayer') else 1
    results[name] = {effect: project(path, layers, magnetic)
                     for effect, magnetic in [('linear_current_driven', False), ('E_times_Bz', True)]}
(here/'allowed.json').write_text(json.dumps(results, indent=2)+'\n')
