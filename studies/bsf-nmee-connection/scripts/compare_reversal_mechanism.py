"""Same-stack spectral controls for the layer NMEE reversal.

H[p,q](k) keeps the zero-field eigenprojectors of endpoint p and assigns the
energy-ordered eigenvalues of q. Physical spin, Zeeman and position operators
are held fixed. The KM velocity is differentiated from H[p,q], rather than
borrowed from either original Hamiltonian. Material electric coupling uses
the unchanged archived Bloch overlaps. These are artificial Hamiltonian
controls, not additive contributions or experimentally realizable structures.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(STUDY), str(STUDY.parents[1] / "papers/jia-2026-nonlin-magnetoelectric-edelstein")]
from bsf_nmee.bilayer_audit import AuditBilayer, AuditSpec, integration_grid
from bsf_nmee.bilayer_protocol import ab_ba_shift
from nmee_lab.clean_core.core import MU_B


def adjoint(a):
    return a.conj().swapaxes(-1, -2)


def transplant(energy, states, target, velocity=None, target_slopes=None):
    """Reject a gauge-dependent assignment inside a degenerate donor group."""
    d = energy[:, None, :] - energy[:, :, None]
    dt = target[:, None, :] - target[:, :, None]
    off = ~np.eye(energy.shape[-1], dtype=bool)
    deg = (abs(d) < 1e-9) & off
    if np.any(abs(dt[deg]) > 1e-7):
        raise ValueError("Target splits a degenerate donor subspace; no unique spectral control")
    h = (states * target[:, None, :]) @ adjoint(states)
    if velocity is None:
        return h
    # A constant occupied/empty spectrum has a projector derivative even at
    # degeneracies within either subspace. Other degenerate derivatives must
    # be obtained from the full Hamiltonian, not individual eigenvectors.
    flat = np.all(target_slopes == 0) and np.all(dt[deg] == 0)
    if deg.any() and not flat:
        raise ValueError("Degenerate velocity block requires a subspace derivative")
    ratio = np.divide(dt, d, out=np.zeros_like(d), where=abs(d) >= 1e-9)
    v = (adjoint(states) @ velocity @ states) * ratio
    ix = np.arange(energy.shape[-1])
    v[:, ix, ix] = target_slopes
    return h, states @ v @ adjoint(states)


def km_response(h, velocity, spins, sz, weights, step, g):
    out = np.zeros(len(spins))
    lower_pairs = np.zeros((4, 4))
    vmax, cmin = -np.inf, np.inf
    lower_k = np.zeros(len(h))
    for sign in (1, -1):
        e, u = np.linalg.eigh(h - sign * step * g * MU_B * sz)
        vmax, cmin = max(vmax, e[:, 3].max()), min(cmin, e[:, 4].min())
        occ, empty = u[:, :, :4], u[:, :, 4:]
        vcv = (adjoint(empty) @ velocity @ occ).swapaxes(1, 2)
        den = (e[:, :4, None] - e[:, None, 4:]) ** 2
        for index, spin in enumerate(spins):
            pairs = -2 * ((adjoint(occ) @ spin @ empty) * vcv).imag / den
            out[index] += sign * (weights @ pairs.sum(axis=(1, 2))) / (2 * step)
            if index == 0:
                lower_pairs += sign * np.einsum('k,kij->ij', weights, pairs) / (2 * step)
                lower_k += sign * pairs.sum(axis=(1, 2)) / (2 * step)
    return out, lower_pairs, vmax, cmin, float(weights @ np.maximum(lower_k, 0)), float(weights @ np.minimum(lower_k, 0))


def run_km(args):
    unit = AuditBilayer(AuditSpec(shift=ab_ba_shift(args.fraction),
                                 staggered_eV=args.mass / 1000, tunneling_eV=1.))
    direction = np.array([.5, np.sqrt(3) / 2])
    spins = [direction[0] * unit.spins[name + '_x'] + direction[1] * unit.spins[name + '_y']
             for name in ['bottom', 'top', 'total']]
    k, weights, _ = integration_grid(unit, *args.grid)
    names = ['weak', 'strong']
    cases = [(p, q) for p in names for q in names] + [(p, 'flat') for p in names]
    rows = {f'{p}_states_{q}_energies': {'C': np.zeros(3), 'pairs': np.zeros((4, 4)),
             'vmax': -np.inf, 'cmin': np.inf, 'positive_k': 0., 'negative_k': 0.} for p, q in cases}
    gap_flat = .030
    models = {name: AuditBilayer(AuditSpec(shift=ab_ba_shift(args.fraction),
                  staggered_eV=args.mass/1000, tunneling_eV=t/1000))
              for name, t in zip(names, args.endpoints)}
    degenerate_errors = []
    degenerate_counts = {f'{p}_{q}': 0 for p, q in cases}

    def control_at(points, donor, receiver):
        e, u = np.linalg.eigh(models[donor].matrices(points))
        target = np.linalg.eigvalsh(models[receiver].matrices(points))
        return transplant(e, u, target)

    def control_velocity(points, donor, receiver, source, target, slopes):
        e, u, v, _ = source
        hc = transplant(e, u, target)
        if donor == receiver:
            return hc, v
        if receiver == 'flat':
            return transplant(e, u, target, v, slopes)
        degenerate = (np.diff(e, axis=1) < 1e-9).any(axis=1)
        vc = np.empty_like(v)
        regular = ~degenerate
        if regular.any():
            _, vc[regular] = transplant(e[regular], u[regular], target[regular],
                                        v[regular], slopes[regular])
        if degenerate.any():
            kp = points[degenerate]
            derivatives = []
            for dk in (1e-5, 5e-6):
                plus = control_at(kp + dk * direction, donor, receiver)
                minus = control_at(kp - dk * direction, donor, receiver)
                derivatives.append((plus - minus) / (2 * dk))
            error = float(np.max(abs(derivatives[0] - derivatives[1])))
            if error > 1e-6:
                raise ValueError(f'Degenerate directional derivative did not converge: {error}')
            degenerate_errors.append(error)
            degenerate_counts[f'{donor}_{receiver}'] += int(degenerate.sum())
            vc[degenerate] = derivatives[1]
        return hc, vc

    min_spacings = []
    for start in range(0, len(k), 512):
        ks, w = k[start:start + 512], weights[start:start + 512]
        hu = unit.matrices(ks)
        vu = sum(direction[j] * unit.matrices(ks, axis=j) for j in range(2))
        sources = {}
        for name, t in zip(names, args.endpoints):
            h, v = hu.copy(), vu.copy()
            for mat in (h, v):
                mat[:, :4, 4:] *= t / 1000
                mat[:, 4:, :4] *= t / 1000
            e, u = np.linalg.eigh(h)
            slopes = np.diagonal(adjoint(u) @ v @ u, axis1=1, axis2=2).real
            sources[name] = (e, u, v, slopes)
            min_spacings.append(float(np.diff(e, axis=1).min()))
        for p, q in cases:
            e, u, v, _ = sources[p]
            if q == 'flat':
                target = np.tile(np.r_[np.full(4, -gap_flat / 2), np.full(4, gap_flat / 2)], (len(e), 1))
                slopes = np.zeros_like(e)
            else:
                target, _, _, slopes = sources[q]
            hc, vc = control_velocity(ks, p, q, sources[p], target, slopes)
            val, pairs, vmax, cmin, pos, neg = km_response(hc, vc, spins, unit.spins['total_z'], w, args.step, unit.mono.g_factor)
            r = rows[f'{p}_states_{q}_energies']
            r['C'] += val
            r['pairs'] += pairs
            r['vmax'], r['cmin'] = max(r['vmax'], vmax), min(r['cmin'], cmin)
            r['positive_k'] += pos
            r['negative_k'] += neg
    for r in rows.values():
        r['common_gap_eV'] = r.pop('cmin') - r.pop('vmax')
        assert r['common_gap_eV'] > 0
        r['layer_sum_error'] = float(r['C'][0] + r['C'][1] - r['C'][2])
        r['C'] = r['C'].tolist()
        r['pairs'] = r['pairs'].tolist()
    # Off-grid derivative validation of the two cross-spectrum Hamiltonians.
    test_k = np.array([[.13, .29], [8.1, -3.2], [-5.8, 12.2]])
    validation_models = [models[name] for name in names]
    errors = []
    for i, j in [(0, 1), (1, 0)]:
        hs = [m.matrices(test_k) for m in validation_models]
        vs = [sum(direction[a] * m.matrices(test_k, axis=a) for a in range(2)) for m in validation_models]
        eig = [np.linalg.eigh(h) for h in hs]
        slopes = np.diagonal(adjoint(eig[j][1]) @ vs[j] @ eig[j][1], axis1=1, axis2=2).real
        _, vc = transplant(*eig[i], eig[j][0], vs[i], slopes)
        displaced = []
        for sign in (1, -1):
            e, u = np.linalg.eigh(validation_models[i].matrices(test_k + sign * 1e-5 * direction))
            eq = np.linalg.eigvalsh(validation_models[j].matrices(test_k + sign * 1e-5 * direction))
            displaced.append(transplant(e, u, eq))
        errors.append(float(np.max(abs(vc - (displaced[0] - displaced[1]) / 2e-5))))
    assert max(errors) < 1e-6
    return {'system': 'KM', 'fraction': args.fraction, 'mass_meV': args.mass,
            'endpoints_t_meV': args.endpoints, 'grid': args.grid, 'nk': len(k), 'step_T': args.step,
            'probe_angle_deg': 60, 'units': '(nm V T)^-1', 'rows': rows,
            'velocity_derivative_errors': errors, 'minimum_band_spacing_eV': min(min_spacings),
            'degenerate_directional_derivative_errors': degenerate_errors,
            'degenerate_directional_derivative_counts': degenerate_counts}


def material_chi(h, spin, z, links, neighbors, cell, nocc, step):
    nk, n = len(h), round(len(h) ** .5)
    values, pair_fields, gaps = [], [], []
    for sign in (1, -1):
        e, u = np.linalg.eigh(h + sign * step * z)
        gaps.append(float(e[:, nocc].min() - e[:, nocc - 1].max()))
        occ, empty = u[:, :, :nocc], u[:, :, nocc:]
        derivative = np.zeros_like(occ)
        for direction in range(4):
            shifted = links[:, direction] @ occ[neighbors[:, direction]]
            overlap = adjoint(occ) @ shifted
            factor = (-1) ** direction * n * cell[direction // 2, 0] / (4 * np.pi)
            derivative += factor * shifted @ np.linalg.inv(overlap)
        svc = adjoint(occ) @ spin @ empty
        rcv = adjoint(empty) @ (1j * derivative)
        pairs = 2 * (svc * rcv.swapaxes(1, 2)).real / (e[:, :nocc, None] - e[:, None, nocc:])
        values.append(float(pairs.sum(axis=(1, 2)).mean()))
        pair_fields.append(pairs)
    c_pairs = (pair_fields[0] - pair_fields[1]) * 6.331927853335898e4 / (2 * step)
    return {'bottom_C_in_1e-9': float((values[0] - values[1]) * 6.331927853335898e4 / (2 * step)),
            'common_gap_eV': min(gaps), 'valence_contributions': c_pairs.mean(axis=0).sum(axis=1).tolist(),
            'k_contributions': c_pairs.sum(axis=(1, 2)).tolist(),
            'positive_k': float(np.maximum(c_pairs.sum(axis=(1, 2)), 0).mean()),
            'negative_k': float(np.minimum(c_pairs.sum(axis=(1, 2)), 0).mean())}


def run_material(args):
    sys.path.insert(0, str(args.data))
    from model import BOTTOM, localized_model
    d = np.load(args.data / 'dp-7269.npz')
    operators = np.load(args.data / 'op-7269.npz')['ops']
    berry = np.load(args.data / 'berry-7269.npz')
    n = round(len(d['k']) ** .5)
    grid, _ = localized_model(d)
    indices = np.rint(d['k'][:, :2] * n).astype(int) % n
    h = np.array([grid[tuple(i)] for i in indices])
    lookup = {tuple(i): j for j, i in enumerate(indices)}
    inverse = [lookup[tuple(-i % n)] for i in indices]
    tr = np.kron(np.eye(22), [[0, 1], [-1, 0]])
    x, z = [(operators[:, i] - tr @ operators[inverse, i].conj() @ tr.T) / 2
            for i in (1, 2)]
    sources = {}
    for name, t in zip(['weak', 'strong'], args.endpoints):
        ht = h * np.where(BOTTOM[:, None] == BOTTOM[None, :], 1., t)
        sources[name] = np.linalg.eigh(ht)
    rows = {}
    for a in sources:
        for b in [*sources, 'flat']:
            e, u = sources[a]
            target = sources[b][0] if b != 'flat' else np.tile(np.r_[np.full(28, -.5), np.full(16, .5)], (len(e), 1))
            ht = transplant(e, u, target)
            err = float(np.max(abs(ht - tr @ ht[inverse].conj() @ tr.T)))
            assert err < 1e-7
            r = material_chi(ht, x, z, berry['mats'], berry['neighbors'], d['cell'], 28, args.step)
            r['TR_residual_eV'] = err
            assert r['common_gap_eV'] > 0
            rows[f'{a}_states_{b}_energies'] = r
    source_hashes = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in
                     [args.data / name for name in ['dp-7269.npz', 'op-7269.npz', 'berry-7269.npz', 'model.py']]}
    return {'system': 'MoSSe 7269 same geometry, 44-state d+p model', 'mesh': n,
            'hopping_endpoints': args.endpoints, 'step_eV': args.step,
            'source_hashes': source_hashes, 'rows': rows,
            'scope': 'Model control, not a new self-consistent DFT calculation or a unique orbital decomposition'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('system', choices=['km', 'mosse'])
    p.add_argument('--data', type=Path)
    p.add_argument('--grid', nargs=3, type=int, default=[41, 45, 72])
    p.add_argument('--fraction', type=float, default=.25)
    p.add_argument('--mass', type=float, default=45.)
    p.add_argument('--endpoints', type=float, nargs=2, required=True)
    p.add_argument('--step', type=float, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if not np.isfinite(a.step) or a.step <= 0:
        p.error('--step must be positive and finite')
    if not np.isfinite(a.endpoints).all() or any(t < 0 for t in a.endpoints):
        p.error('--endpoints must be finite and nonnegative')
    if a.system == 'mosse' and a.data is None:
        p.error('--data is required for the material model')
    result = run_km(a) if a.system == 'km' else run_material(a)
    result['method'] = 'zero-field spectral transplantation, energy-order correspondence; add physical Zeeman perturbation afterwards'
    result['interpretation'] = ('Finite Hamiltonian controls, not an additive energy/geometry decomposition. '
        'Energy order fixes the correspondence. Band-resolved diagnostics within degenerate '
        'subspaces are not unique; only complete subspace sums are physical. '
        'Flat controls alter all transition energies and are not predictions for a realizable band structure.')
    for name, row in result['rows'].items():
        if name.endswith('_flat_energies'):
            row.pop('pairs', None)
            row.pop('valence_contributions', None)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    for name, row in result['rows'].items():
        print(name, row.get('C', row.get('bottom_C_in_1e-9')), 'gap', row['common_gap_eV'], flush=True)


if __name__ == '__main__':
    main()
