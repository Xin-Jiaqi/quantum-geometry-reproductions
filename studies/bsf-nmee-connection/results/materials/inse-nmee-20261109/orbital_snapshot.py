"""Gamma band-pair characters in the existing orthonormalized In-s,p / Se-p trial basis.

Descriptive orbital character, not an additive decomposition of NMEE and not
an orbital-selective causal intervention. Kramers partners are summed together.
"""
from pathlib import Path
import os, socket, json
import numpy as np
from ase.io import read
from analyze import embedding
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super', 'cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root = Path.cwd(); d = np.load(root/'A0/A/cache.npz'); at = read(root/'A0/A/seed/POSCAR')
labels = []; lower = []
for atom in at:
    # The existing real-harmonic trial order is s then p(m=-1,0,+1), spin fastest.
    groups = ['In_s', 'In_pxy', 'In_pz', 'In_pxy'] if atom.symbol == 'In' else ['Se_pxy', 'Se_pz', 'Se_pxy']
    for name in groups:
        labels.extend([name]*2); lower.extend([atom.position[2] < at.cell[2, 2]/2]*2)
labels = np.array(labels); lower = np.array(lower)
assert len(labels) == len(d['orb_bottom']) and np.array_equal(lower, d['orb_bottom'])
ig = int(np.argmin(np.linalg.norm(d['k'], axis=1))); assert np.max(abs(d['k'][ig])) < 1e-8
no = int(d['nocc']); vbm = d['energies'][:, no-1].max()
result = {'scope': __doc__, 'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(),
          'energy_reference': 'global sampled VBM at zero gate', 'rows': []}
for nb in (128, 144, 160):
    for frozen in (4, 8, 12):
        lift, sv = embedding(d, frozen, nb)
        for first in range(no-8, no+4, 2):
            weights = abs(lift[ig, first:first+2])**2
            capture = weights.sum(axis=1)
            assert np.max(abs(capture-1)) < 1e-8
            w = weights.mean(axis=0)
            result['rows'].append({'nb': nb, 'frozen_empty': frozen, 'bands_1_based': [first+1, first+2],
                'energy_eV': (d['energies'][ig, first:first+2]-vbm).tolist(),
                'character': {name: float(w[labels == name].sum()) for name in sorted(set(labels))},
                'lower_weight': float(w[lower].sum()), 'min_trial_singular': sv})
(root/'orbital_snapshot.json').write_text(json.dumps(result, indent=2))
print('ORBITAL_SNAPSHOT_COMPLETE', flush=True)
