"""Spatial selection rules for the actual parent monolayer; no DFT or switching model."""
from pathlib import Path
import json
import numpy as np
import spglib
from ase.io import read

p = Path(__file__).resolve().parent
a = read(p.parent/'inse-nmee-20261109/POSCAR.mono')
sym = spglib.get_symmetry((a.cell.array, a.get_scaled_positions(), a.numbers), symprec=1e-5)
L = a.cell.array.T
ops = [L@r@np.linalg.inv(L) for r in sym['rotations']]
# Spin and B are axial: det(R)^2=1. The spatial rank-three representation
# consequently equals the polar one, without equating microscopic mechanisms.
matrix = np.vstack([np.kron(np.kron(r,r),r)-np.eye(27) for r in ops])
_, sv, vh = np.linalg.svd(matrix); invariant = vh[sum(sv>1e-8):].T
allowed = ['xyz'[i]+'xyz'[j]+'xyz'[k] for i in range(3) for j in range(3) for k in range(3)
           if np.linalg.norm(invariant[9*i+3*j+k])>1e-8]
assert 'xxz' not in allowed and invariant.shape[1] == 1
out = {'operation_count': len(ops), 'tensor': 's_i=Gamma_ijk E_j B_k; two axial factors cancel det(R)^2 under spatial operations',
       'invariant_dimension': invariant.shape[1], 'allowed_indices_in_actual_cartesian_frame': allowed,
       'xxz_allowed': False,
       'layer_force_argument': 'C3 leaves F_z invariant; no nonzero invariant in-plane layer force vector. This is not a positive-curvature or switching-threshold test.',
       'FQFE_scope': 'An isolated D3h endpoint keeps xxz forbidden at zero gate. Allowed all-in-plane Gamma components reverse under an inversion-related endpoint. A vertical gate cannot be inferred to select these endpoints; FQFE changes internal sublattice registry, not interlayer displacement.'}
(p/'symmetry_checks.json').write_text(json.dumps(out,indent=2))
