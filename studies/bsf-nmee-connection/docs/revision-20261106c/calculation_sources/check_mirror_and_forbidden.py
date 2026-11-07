"""Check the crystallographic mirror and forbidden response near the first zero."""
from runtime import STUDY, OUTPUT
import sys
import json
from types import SimpleNamespace
import numpy as np

sys.path[:0] = [str(STUDY/'scripts'), str(STUDY),
               str(STUDY.parents[1]/'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from run_fixed_stack_layer_cxx import response
from bsf_nmee.bilayer_audit import AuditBilayer, AuditSpec, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift

normal = np.array([.5, np.sqrt(3)/2])
transverse = np.array([-np.sqrt(3)/2, .5])
reflection = np.eye(2) - 2*np.outer(normal, normal)
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
spin_mirror = np.kron(np.eye(4), normal[0]*sigma_x + normal[1]*sigma_y)
k_test = np.array([[.21, .33], [1.17, -.49], [-.73, 1.21]])
mirror_rows = []
for fraction in [0, .2, .25, .3, 1]:
    model = AuditBilayer(AuditSpec(shift=ab_ba_shift(fraction), staggered_eV=.045,
                                  tunneling_eV=.010))
    error = np.max(abs(model.matrices(k_test@reflection.T) -
                       spin_mirror@model.matrices(k_test)@spin_mirror.conj().T))
    assert error < 1e-12, error
    mirror_rows.append({'lambda': fraction, 'max_error_eV': float(error)})

model = AuditBilayer(AuditSpec(shift=ab_ba_shift(.25), staggered_eV=.045,
                              tunneling_eV=1.))
k, weights, _ = integration_grid(model, 41, 45, 72)
h = np.concatenate([model.matrices(k[i:i+512]) for i in range(0, len(k), 512)])
v = np.concatenate([normal[0]*model.matrices(k[i:i+512], axis=0) +
                    normal[1]*model.matrices(k[i:i+512], axis=1)
                    for i in range(0, len(k), 512)])
result = {'mirror': mirror_rows, 'lambda': .25, 'mass_meV': 45,
          'grid': [41, 45, 72], 'step_T': .0625, 't_meV': 10,
          'U_meV': 4.75, 'probe_direction_native_XY': normal.tolist(), 'response': {}}
for label, vector in [('longitudinal', normal), ('transverse', transverse), ('out_of_plane', None)]:
    spins = {}
    for layer in ['bottom', 'top', 'total']:
        spins[layer+'_x'] = (model.spins[layer+'_z'] if vector is None else
                            vector[0]*model.spins[layer+'_x'] + vector[1]*model.spins[layer+'_y'])
    spins['total_z'] = model.spins['total_z']  # Global Zeeman vertex is never rotated out of z.
    probe = SimpleNamespace(mono=model.mono, spins=spins)
    row = response(probe, h, v, weights, 10, 4.75, .0625)
    assert row['common_gap_meV'] > 0
    assert abs(row['layer_sum_error']) < 1e-12
    if label != 'longitudinal':
        assert max(abs(row[x]) for x in ['bottom', 'top', 'total']) < 1e-10
    result['response'][label] = row
    print(label, row['bottom'], row['top'], row['total'], flush=True)
(OUTPUT/'forbidden_at_compensation.json').write_text(json.dumps(result, indent=2)+'\n')
