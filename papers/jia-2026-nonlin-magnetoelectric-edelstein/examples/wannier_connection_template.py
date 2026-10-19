"""TEMPLATE, not a claim of completed ab initio NMEE.

Required USER data (not supplied in this archive):
  - physical 2D in-plane lattice vectors in nm;
  - spin_onsite.npz matrices 'x', 'y', 'z' in Pauli-sigma convention;
  - actual orbital-layer mapping for the *same* Wannier ordering;
  - position-operator/Berry-connection information if quantitative physics needed.

This code opts into zero connection as an EXPLORATORY test ONLY.
"""
from pathlib import Path
import numpy as np
from nmee_lab import read_wannier90_hr,to_clean_model,mesh,evaluate_clean,CleanSpec

# For a real material replace the placeholder paths, lattice vectors and layer map.
# lattice_nm = np.array([[..., ...], [..., ...]], dtype=float)
# tb = read_wannier90_hr('wannier90_hr.dat', lattice_nm)
# with np.load('spin_onsite.npz') as f:
#     pauli_spin = {direction: f[direction] for direction in ('x','y','z')}
# layers = {'top': [...], 'bottom': [...]} # 0-based orbital indices
# model = to_clean_model(tb,spin_onsite=pauli_spin,layer_indices=layers,
#                        acknowledge_zero_connection=True)
# k=mesh(model.reciprocal,81)
# for a in ('top_x','bottom_x','x'):
#     print(a,evaluate_clean(model,k,CleanSpec(alpha=a,gamma='z',beta=0,mu_eV=0.,temperature_K=10.)))
# Note: do not interpret these as calibrated predictions without the connections.
