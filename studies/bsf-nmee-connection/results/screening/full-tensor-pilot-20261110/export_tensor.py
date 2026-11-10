"""Full-KS PAW spin vertices and 2D links. IMR science requires PBS compute node.

No Wannier interpolation/downfolding; MMN supplies neighboring full-KS overlaps.
Both spatial layers are integrated independently, including PAW augmentation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import socket
import sys
import time
import numpy as np


def adj(x):
    return x.conj().swapaxes(-1, -2)


def spin_blocks(w, projected):
    a, b = w
    c, d = projected
    return np.array([adj(a) @ d + adj(b) @ c,
                     -1j * adj(a) @ d + 1j * adj(b) @ c,
                     adj(a) @ c - adj(b) @ d])


def masked_wave(w, g, mask):
    fourier = np.fft.ifft(mask)
    result = np.empty_like(w)
    _, groups = np.unique(g[:, :2], axis=0, return_inverse=True)
    for group in np.unique(groups):
        ids = np.flatnonzero(groups == group)
        gz = g[ids, 2]
        block = fourier[(gz[None, :] - gz[:, None]) % len(mask)]
        for s in (0, 1):
            result[s, ids] = block @ w[s, ids]
    return result


def vertices(w, g, pro, masks, atom_lower):
    total = spin_blocks(w, w)
    layers = [spin_blocks(w, masked_wave(w, g, m)) for m in masks]
    metric = sum(adj(v) @ v for v in w)
    for ia, typ in enumerate(pro.element_idx):
        p = pro.qproj[typ] * pro.cqfak[typ][:, None] * pro.crexp[None, :, ia]
        beta = np.array([p @ v for v in w])
        qb = np.array([pro.pawpp[typ].get_Qij() @ v for v in beta])
        aug = spin_blocks(beta, qb)
        total += aug
        layers[0 if atom_lower[ia] else 1] += aug
        metric += sum(adj(v) @ qv for v, qv in zip(beta, qb))
    all_spin = np.array([total, *layers])
    herm_error = float(np.max(np.abs(all_spin - adj(all_spin))))
    return (all_spin + adj(all_spin)) / 2, metric, herm_error


def export(seed, output, layer_mode, mesh, cut_override=None):
    from ase.io import read
    from paw import pawpotcar, nonlq
    from vaspwfc import vaspwfc
    from spin_transitions import fft_z, spin_and_metric
    start = time.monotonic()
    atoms = read(seed / 'POSCAR')
    if np.max(abs(atoms.cell.array[:2, 2])) > 1e-8 or np.max(abs(atoms.cell.array[2, :2])) > 1e-8:
        raise ValueError('Slab must lie in Cartesian xy with orthogonal vacuum z')
    wave = vaspwfc(str(seed / 'WAVECAR'), lsorbit=True)
    pots = [pawpotcar(s) for s in (seed / 'POTCAR').read_text().split('End of Dataset')[:-1]]
    nb = wave._nbands
    k = wave._kvecs
    nk = len(k)
    if nk != int(np.prod(mesh)):
        raise ValueError('Full unreduced N1*N2 k mesh required')
    nocc = int((seed / 'EIGENVAL').read_text().splitlines()[5].split()[0])
    z = atoms.get_scaled_positions()[:, 2]
    cut = .5
    if layer_mode == 'bilayer':
        zs = np.sort(z)
        if len(zs) % 2:
            raise ValueError('Homobilayer needs even atom count')
        cut = float(np.mean(zs[len(zs)//2-1:len(zs)//2+1]))
        if zs[len(zs)//2] - zs[len(zs)//2-1] <= 0:
            raise ValueError('No interlayer gap')
    if cut_override is not None:
        cut = float(cut_override)
    if not np.isfinite(cut):
        raise ValueError('Layer cut must be finite')
    cut %= 1
    if layer_mode == 'bilayer' and not zs[len(zs)//2-1] < cut < zs[len(zs)//2]:
        raise ValueError('Layer cut must be inside central interlayer gap')
    if layer_mode == 'bilayer' and (np.any(np.isclose((z-cut) % 1, 0, atol=1e-10, rtol=0)) or np.any(np.isclose((z-cut) % 1, .5, atol=1e-10, rtol=0))):
        raise ValueError('Atom intersects spatial layer boundary')
    nz = fft_z(seed / 'CHGCAR')
    t = (np.arange(nz) / nz - cut) % 1
    lower = (t >= .5).astype(float)
    lower[np.isclose(t, 0, atol=1e-12, rtol=0) | np.isclose(t, .5, atol=1e-12, rtol=0)] = .5
    masks = [lower, 1-lower]
    atom_lower = ((z-cut) % 1 >= .5)
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    links = np.zeros((nk, 4, nb, nb), complex)
    neighbors = np.full((nk, 4), -1, int)
    with (seed / 'wannier90.mmn').open() as f:
        next(f)
        bands, nn, nh = map(int, next(f).split())
        if (bands, nn) != (nb, nk):
            raise ValueError('MMN band/grid mismatch')
        for _ in range(nk * nh):
            i, j, *shift = map(int, next(f).split())
            i -= 1
            j -= 1
            block = np.fromstring(' '.join(next(f) for _ in range(nb*nb)), sep=' ').reshape(-1, 2)
            vector = k[j] + shift - k[i]
            scaled = vector[:2] * mesh
            b = tuple(np.rint(scaled).astype(int))
            if abs(vector[2]) < 1e-7 and np.max(abs(scaled-np.rint(scaled))) < 1e-7 and b in steps:
                q = steps.index(b)
                if neighbors[i, q] != -1:
                    raise ValueError('Duplicate directional link')
                neighbors[i, q] = j
                links[i, q] = (block[:, 0]+1j*block[:, 1]).reshape(nb, nb).T
    if np.any(neighbors < 0):
        raise ValueError('Missing directional full-KS link')
    reverse = max(float(np.max(abs(links[i, q]-adj(links[j, q ^ 1]))))
                  for i in range(nk) for q, j in enumerate(neighbors[i]))
    if reverse > 1e-6:
        raise ValueError('MMN reverse-link mismatch')
    # Independent PAW overlap check: WAVECAR and MMN must share a gauge.
    i, direction = 0, 0
    j = neighbors[i, direction]
    shift = np.rint(k[i] + [1/mesh[0], 0, 0] - k[j]).astype(int)
    gs, ws, ps = [], [], []
    for ix in (i, j):
        gs.append(wave.gvectors(ix+1))
        ws.append(np.stack([wave.readBandCoeff(ikpt=ix+1, iband=b+1).reshape(2, -1) for b in range(nb)], axis=-1))
        ps.append(nonlq(atoms, wave._encut, potcar=pots, k=k[ix]))
    jmap = {tuple(g): q for q, g in enumerate(gs[1])}
    ii = [q for q, g in enumerate(gs[0]) if tuple(g+shift) in jmap]
    jj = [jmap[tuple(gs[0][q]+shift)] for q in ii]
    overlap = sum(adj(ws[0][s, ii]) @ ws[1][s, jj] for s in (0, 1))
    for ia, typ in enumerate(ps[0].element_idx):
        proj = [v.qproj[typ]*v.cqfak[typ][:, None]*v.crexp[None, :, ia] for v in ps]
        phase = np.exp(-2j*np.pi*shift @ atoms.get_scaled_positions()[ia])
        for spin in (0, 1):
            overlap += phase * adj(proj[0] @ ws[0][spin]) @ pots[typ].get_Qij() @ (proj[1] @ ws[1][spin])
    window = slice(max(0, nocc-8), min(nb, nocc+8))
    gauge_error = float(np.linalg.norm((overlap-links[i, direction])[window, window]) / np.linalg.norm(links[i, direction][window, window]))
    if gauge_error > .08:
        raise ValueError(f'Wave/link gauge mismatch: {gauge_error}')
    spins, metrics, herm, old_errors = [], [], [], []
    for i, kv in enumerate(k):
        g = wave.gvectors(i+1)
        w = np.stack([wave.readBandCoeff(ikpt=i+1, iband=b+1).reshape(2, -1) for b in range(nb)], axis=-1)
        pro = nonlq(atoms, wave._encut, potcar=pots, k=kv)
        if not np.array_equal(g, pro.Gvec):
            raise ValueError('PAW plane-wave ordering mismatch')
        s, metric, he = vertices(w, g, pro, masks, atom_lower)
        err = float(np.max(abs(metric-np.eye(nb))))
        if err > 1e-4 or he > 1e-7:
            raise ValueError(f'PAW metric/hermiticity failed: {err}, {he}')
        # Old implementation independently confirms the x-spin lower vertex.
        if i == 0:
            _, legacy = spin_and_metric(w, pro, g, masks[0], atom_lower)
            old_errors.append(float(np.max(abs(legacy-s[1, 0]))))
        spins.append(s)
        metrics.append(err)
        herm.append(he)
        print('EXPORT_TENSOR', i+1, nk, flush=True)
    spins = np.array(spins)
    summed = float(np.max(abs(spins[:, 0]-spins[:, 1]-spins[:, 2])))
    if summed > 1e-7 or max(old_errors) > 1e-9:
        raise ValueError('Independent layer integration or legacy sx mismatch')
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'tensor_cache.npz', energies=wave._bands[0],
                        Z=spins[:, 0, 2], spin_tensor=spins, links=links, neighbors=neighbors,
                        mesh=mesh, k=k, cell=atoms.cell.array, nocc=nocc)
    record = dict(nk=nk, nb=nb, nocc=nocc, mesh=mesh.tolist(), layer_mode=layer_mode,
                  cut_fractional_z=cut, layer_semantics='physical layers' if layer_mode=='bilayer' else 'only total is physical isolated-layer response',
                  metric_max=max(metrics), hermiticity_max=max(herm), layer_sum_operator_error=summed,
                  legacy_lower_sx_operator_error=max(old_errors), link_reverse_error=reverse, wave_link_gauge_error=gauge_error,
                  export_seconds=time.monotonic()-start, node=socket.gethostname(), job=os.environ['PBS_JOBID'],
                  source_seed=str(seed), input_sha256={x:hashlib.sha256((seed/x).read_bytes()).hexdigest() for x in ['POSCAR','INCAR','KPOINTS']})
    (output / 'export.json').write_text(json.dumps(record, indent=2))
    return record


if __name__ == '__main__':
    if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super', 'cms-ssh')):
        raise SystemExit('PBS compute allocation required')
    p = argparse.ArgumentParser()
    p.add_argument('seed', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--layer-mode', choices=['bilayer','mono'], required=True)
    p.add_argument('--mesh', nargs=2, type=int, required=True)
    p.add_argument('--cut', type=float, default=None)
    args = p.parse_args()
    parent = Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
    sys.path[:0] = [str(parent), str(parent/'pawlib')]
    export(args.seed, args.output, args.layer_mode, np.array(args.mesh), args.cut)
