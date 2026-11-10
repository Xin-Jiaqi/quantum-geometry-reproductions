"""Local structural/representation audit; no DFT or response integration."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
import numpy as np
import spglib
from ase.io import read

HERE = Path(__file__).resolve().parent
PLAN = HERE.parent / 'bidb-pilot-plan-20261110'
TOOLKIT = Path('/Users/jiaqi/Desktop/worksofJiaqi/project-x/tmp/repos/group-theory-operations-toolkit-nmee-v018')
sys.path.insert(0, str(TOOLKIT / 'src'))
from group_theory_operations import equivariant_map_basis
from group_theory_operations.stacking import partition_left_cosets, equivalent_interface_translation


def symmetry(atoms):
    return spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-3)


def tensor_check(atoms, layered):
    ds = symmetry(atoms)
    lattice = atoms.cell.array.T
    output, inputs = [], []
    for rf in ds.rotations:
        r = lattice @ rf @ np.linalg.inv(lattice)
        axial = np.linalg.det(r) * r
        permutation = np.eye(2) if r[2, 2] > 0 else np.array([[0, 1], [1, 0]])
        output.append(np.kron(permutation, axial) if layered else axial)
        inputs.append(r[:2, :2] * axial[2, 2])
    checks = {}
    for gate in (False, True):
        ins = [v * (rf[2, 2] if gate else 1) for v, rf in zip(inputs, ds.rotations)]
        raw = equivariant_map_basis([v.tolist() for v in output], [v.tolist() for v in ins])
        basis = np.array(raw).reshape(-1, 2 if layered else 1, 3, 2)
        for b in basis:
            assert np.max(np.abs(b[:, 0, 0] - b[:, 1, 1])) < 1e-8
            assert np.max(np.abs(b[:, [0, 1, 2, 2], [1, 0, 0, 1]])) < 1e-8
            if ds.pointgroup == '-3m':
                assert np.max(np.abs(b[0] - b[1] if gate else b[0] + b[1])) < 1e-8
        checks['gate_derivative' if gate else 'zero_gate'] = {'dimension': len(basis), 'basis': basis.tolist()}
    return checks


inventory = json.loads((PLAN / 'inventory.json').read_text())
report = {
    'toolkit_commit': inventory['toolkit_commit'],
    'coordinates': 'a1=a(1,0,0), a2=a(-1/2,sqrt(3)/2,0); +z lower to upper; lower-layer mean z reset to0; in-plane input origin unchanged',
    'seitz_convention': 'r_upper=O r_lower+tau1*a1+tau2*a2+d*z_hat; d is layer-centroid separation, not interface gap; no modulo reduction along nonperiodic z',
    'representative_caution': 'Selected O is a representative of O G_S; not a unique stacking name. Actual bilayer symmetry must map both layers, not merely lower onto upper.',
    'material_geometry_scope': 'BiDB reference structures at symprec1e-3A. A relaxed structure needs its own symmetry and displacement check; do not assert exact same-mother mapping after independent intralayer relaxation.',
    'monolayers': {}, 'bilayers': [],
}
for mid in (7265, 1929, 862):
    atoms = read(PLAN / f'structures/mono-{mid}.vasp')
    ds = symmetry(atoms)
    parent = spglib.get_symmetry((atoms.cell.array, [[0, 0, .5]], [1]), symprec=1e-3)
    cosets = partition_left_cosets(parent['rotations'], ds.rotations)
    mirror = np.diag([1, 1, -1])
    report['monolayers'][str(mid)] = {
        'pointgroup': ds.pointgroup, 'orientation_left_coset_count': len(cosets),
        'm001_in_monolayer_pointgroup': any(np.array_equal(r, mirror) for r in ds.rotations),
        'tensor': tensor_check(atoms, False),
    }
for row in inventory['bilayers']:
    if row['id'] not in (7266, 7268, 7269, 7270, 1930, 1931, 864, 863):
        continue
    path = PLAN / f"structures/bilayer-{row['id']}.vasp"
    atoms = read(path)
    ds = symmetry(atoms)
    f = atoms.get_scaled_positions(wrap=False)
    order = np.argsort(atoms.positions[:, 2])
    lower, upper = np.array_split(order, 2)
    f[:, 2] -= f[lower, 2].mean()
    conn = row['connection']
    rotation = np.array(conn['fractional_rotation'])
    tau = [float(Fraction(x).limit_denominator(12)) for x in conn['inplane_translation_mod1']]
    d = float((f[upper, 2].mean() - f[lower, 2].mean()) * atoms.cell[2, 2])
    translation = np.array([*tau, d / atoms.cell[2, 2]])
    errors, matches = [], []
    for i in lower:
        delta = f[upper] - (rotation @ f[i] + translation)
        delta[:, :2] -= np.rint(delta[:, :2])
        distance = np.linalg.norm(delta @ atoms.cell.array, axis=1)
        distance[atoms.numbers[upper] != atoms.numbers[i]] = np.inf
        j = int(distance.argmin())
        errors.append(float(distance[j]))
        matches.append(j)
    assert max(errors) < 1e-6 and len(set(matches)) == len(lower)
    report['bilayers'].append({
        'id': row['id'], 'material': row['parent_formula'], 'parent_id': row['parent_id'],
        'pointgroup': ds.pointgroup, 'stacking_point_representative': conn['representative'],
        'tau_fractional': [str(Fraction(x).limit_denominator(12)) for x in tau],
        'd_centroid_A': d, 'lower_to_upper_mapping_error_A': max(errors),
        'layer_preserving_operations': int(sum(r[2, 2] > 0 for r in ds.rotations)),
        'layer_exchanging_operations': int(sum(r[2, 2] < 0 for r in ds.rotations)),
        'tensor': tensor_check(atoms, True), 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
    })
report['mirror_partner_tau_A_to_B'] = equivalent_interface_translation([1/3, 2/3], np.diag([1, 1, -1]))
(HERE / 'structural_tensor_audit.json').write_text(json.dumps(report, indent=2) + '\n')
print('Verified eight Seitz mappings, parent left cosets and layer-resolved EB/gate tensor constraints.')

# KM has a separate lattice embedding: its stored displacement moves the lower layer.
repo = next(p for p in HERE.parents if (p / 'studies/bsf-nmee-connection/bsf_nmee').is_dir())
sys.path[:0] = [str(repo / 'studies/bsf-nmee-connection'), str(repo / 'papers/jia-2026-nonlin-magnetoelectric-edelstein')]
from bsf_nmee.bilayer_audit import AuditSpec, AuditBilayer
from bsf_nmee.bilayer_protocol import ab_ba_shift
from nmee_lab.clean_core.models import SX, SY
n = np.array([.5, np.sqrt(3)/2])
reflection = np.eye(2) - 2*np.outer(n, n)
unitary = np.kron(np.eye(4), n[0]*SX+n[1]*SY)
kpoints = np.array([[.13, .29], [8.1, -3.2], [-5.8, 12.2]])
km_rows = []
for stacking in ('direct', 'mirror'):
    model = AuditBilayer(AuditSpec(stacking=stacking, shift=ab_ba_shift(.25), staggered_eV=.045, tunneling_eV=.010))
    tau = -np.array(ab_ba_shift(.25)) @ model.a
    counts = {1: 0, -1: 0}
    for j in range(3):
        angle = 2*np.pi*j/3
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        for h in (rotation, rotation @ np.diag([-1, 1])):
            for sign in (1, -1):
                if sign == -1 and stacking == 'direct':
                    continue  # identical C3v monolayers have no normal-reversing symmetry
                mismatch = (h @ tau-sign*tau) @ np.linalg.inv(model.a)
                counts[sign] += int(np.max(np.abs(mismatch-np.rint(mismatch))) < 1e-8)
    residual = float(np.max(np.abs(model.matrices(kpoints @ reflection.T)-unitary @ model.matrices(kpoints) @ unitary.conj().T)))
    assert counts == {1: 2, -1: 0} and residual < 1e-9
    km_rows.append({'stacking': stacking, 'pointgroup': 'Cs(m)', 'lambda': .25,
                    'layer_preserving_count': counts[1], 'layer_exchanging_count': counts[-1],
                    'vertical_mirror_residual_eV': residual,
                    'lower_to_upper_tau_code_60degree_basis': (-np.array(ab_ba_shift(.25))).tolist()})
(HERE / 'km_coordinate_check.json').write_text(json.dumps({
    'scope': 'BSF registry conditions and three generic k-point Hamiltonian mirror checks, no response integration',
    'native_basis_rows_nm': model.a.tolist(), 'probe_xprime_native_xy': n.tolist(),
    'point_representatives_lower_to_upper': {'direct': 'E', 'mirror': 'm001'},
    'source_sha256': {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [repo/'studies/bsf-nmee-connection/bsf_nmee/bilayer_audit.py', repo/'studies/bsf-nmee-connection/bsf_nmee/bilayer_protocol.py', repo/'studies/bsf-nmee-connection/bsf_nmee/sliding_tight_binding.py']},
    'rows': km_rows}, indent=2)+'\n')
