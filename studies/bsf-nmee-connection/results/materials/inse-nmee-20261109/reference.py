"""Matched InSe smooth planar electrostatic reference at three screened gate fields.

The potential comes from the actual self-consistent bilayer minus its matched
isolated layer at the SAME applied field. This is a specified frozen-potential
reference, not a unique or additive partition of interface physics. No new SCF.
PAW augmentation uses the locally constant smooth potential approximation.
"""
from pathlib import Path
import json, os, socket, sys
import numpy as np
from ase.io import read
from analyze import Response

if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super', 'cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root = Path.cwd()
parent = Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0] = [str(parent), str(parent/'pawlib')]
from paw import pawpotcar, nonlq
from vaspwfc import vaspwfc

def plane(path):
    with path.open() as f:
        h = [next(f) for _ in range(8)]
        for _ in range(sum(map(int, h[6].split()))):
            next(f)
        nx, ny, nz = map(int, next(s for s in f if s.strip()).split())
        vals = []
        while len(vals) < nx*ny*nz:
            vals.extend(map(float, next(f).split()))
    assert len(vals) == nx*ny*nz
    return np.array(vals).reshape(nz, ny, nx).mean(axis=(1, 2))

def shift(v, frac):
    return np.fft.ifft(np.fft.fft(v)*np.exp(2j*np.pi*np.fft.fftfreq(len(v))*len(v)*frac)).real

def smooth(v, sigma, length):
    return np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(len(v))*len(v)/length*sigma)**2)).real

planar = {}
result = {'job': os.environ['PBS_JOBID'], 'node': socket.gethostname(), 'scope': __doc__,
          'units': '1e-9 (nm V T)^-1, Pauli spin density', 'matched_intralayer_geometry': True,
          'smooth_sigma_A': .2, 'points': {}, 'complete': False}
for bilabel, monolabel, field in [('A0', 'mono0', 0.), ('Am', 'monom', -.1), ('Ap', 'monop', .1)]:
    base = root/monolabel/'mono'/'seed'
    ab = read(root/bilabel/'A'/'seed'/'POSCAR')
    at = read(base/'POSCAR')
    assert np.max(abs(ab.cell.array-at.cell.array)) < 1e-8
    vm = plane(base/'LOCPOT'); vb = plane(root/bilabel/'A'/'seed'/'LOCPOT')
    assert len(vm) == len(vb)
    length = at.cell[2, 2]; nz = len(vm); frac = np.arange(nz)/nz; z = (frac-.5)*length
    planar[bilabel] = {'F': field, 'z_A': (frac*length).tolist(),
                       'bilayer_potential_eV': vb.tolist(), 'mono_potential_eV': vm.tolist()}
    (root/'planar.json').write_text(json.dumps(planar, indent=2))
    profiles = []; specs = []
    # Match both layers by a rigid translation; no reflection of the spin basis.
    for name, mask in [('lower', ab.positions[:, 2] < length/2), ('upper', ab.positions[:, 2] > length/2)]:
        layer = ab[mask]
        center = layer.positions[:, 2].mean()
        shifted = layer.positions.copy(); shifted[:, 2] += length/2-center
        shiftxy = at.positions[0, :2]-shifted[0, :2]
        shifted[:, :2] += shiftxy
        delta = (shifted-at.positions)@np.linalg.inv(at.cell.array)
        delta[:, :2] -= np.rint(delta[:, :2])
        assert layer.get_chemical_symbols() == at.get_chemical_symbols()
        assert np.max(abs(delta@at.cell.array)) < 1e-6
        halfwidth = np.max(abs(layer.positions[:, 2]-center))
        diff = shift(vb, center/length-.5)-vm
        diff -= diff[abs(z) < halfwidth].mean()
        diff = smooth(diff, .20, length)
        start = halfwidth+.4
        for margin in (1.2, 1.8):
            end = halfwidth+margin
            u = np.clip((abs(z)-start)/(end-start), 0, 1)
            profiles.append(diff*.5*(1+np.cos(np.pi*u)))
            specs.append(dict(layer=name, taper_start_A=float(start), taper_end_A=float(end)))
    profiles = np.array(profiles); vfft = np.fft.ifft(profiles, axis=1)
    wave = vaspwfc(str(base/'WAVECAR'), lsorbit=True); nb = wave._nbands
    pots = [pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]]
    operators = []
    for ik, k in enumerate(wave._kvecs):
        g = wave.gvectors(ik+1)
        w = np.stack([wave.readBandCoeff(ikpt=ik+1, iband=b+1).reshape(2, -1) for b in range(nb)], axis=-1)
        pr = nonlq(at, wave._encut, potcar=pots, k=k)
        out = np.zeros((len(profiles), nb, nb), complex)
        _, groups = np.unique(g[:, :2], axis=0, return_inverse=True)
        for group in np.unique(groups):
            ids = np.where(groups == group)[0]; gz = g[ids, 2]
            block = vfft[:, (gz[None, :]-gz[:, None]) % nz]
            for spin in (0, 1):
                out += w[spin, ids].conj().T@(block@w[spin, ids])
        for ia, typ in enumerate(pr.element_idx):
            op = pr.qproj[typ]*pr.cqfak[typ][:, None]*pr.crexp[None, :, ia]
            betas = [op@ww for ww in w]; q = pots[typ].get_Qij()
            values = np.array([np.interp(at.get_scaled_positions()[ia, 2], frac, v) for v in profiles])
            out += values[:, None, None]*sum(b.conj().T@q@b for b in betas)
        operators.append((out+out.conj().transpose(0, 2, 1))/2)
        print('REFERENCE_OPERATOR', bilabel, ik+1, len(wave._kvecs), flush=True)
    operators = np.array(operators).transpose(1, 0, 2, 3)
    np.savez_compressed(root/f'reference_{bilabel}.npz', operators=operators, profiles=profiles, z_A=z)
    d = np.load(root/monolabel/'mono'/'cache.npz'); rows = []
    for cutoff in (64, 72, 80):
        r = Response(d, cutoff)
        for i, spec in enumerate(specs):
            ans = r.response(r.h+operators[i, :, :cutoff, :cutoff])
            rows.append(dict(nb=cutoff, **spec, C_layer=ans['C_total_lower_upper'][0], **ans))
    result['points'][bilabel] = dict(F_external_V_nm=field, rows=rows)
    (root/'reference.json').write_text(json.dumps(result, indent=2))
gate_checks = {}
for system in ('bilayer', 'mono'):
    dv = (np.array(planar['Ap'][system+'_potential_eV'])-planar['Am'][system+'_potential_eV'])/.2
    zz = np.array(planar['Ap']['z_A']); vacuum = (zz > 2.) & (zz < 4.)
    fit = np.polyfit(zz[vacuum], dv[vacuum], 1)
    gate_checks[system] = {'vacuum_dV_dz_per_F': float(fit[0]),
                           'vacuum_fit_residual': float(np.std(dv[vacuum]-np.polyval(fit, zz[vacuum])))}
    if system == 'bilayer':
        # Average near each of the two In planes in each layer. This is a local
        # lever arm diagnostic, not a unique dielectric constant.
        planes = sorted(p[2] for s, p in zip(ab.get_chemical_symbols(), ab.positions) if s == 'In')
        values = [float(dv[abs(zz-zp) < .25].mean()) for zp in planes]
        gate_checks[system]['In_planes_z_A'] = planes
        gate_checks[system]['In_plane_dV_per_F'] = values
        gate_checks[system]['interlayer_lever_eV_per_V_nm'] = float(np.mean(values[2:])-np.mean(values[:2]))
        gate_checks[system]['bare_linear_lever_eV_per_V_nm'] = float((np.mean(planes[2:])-np.mean(planes[:2]))/10)
result['gate_potential_checks'] = gate_checks
result['complete'] = True
(root/'reference.json').write_text(json.dumps(result, indent=2))
print('REFERENCE_DONE', flush=True)
