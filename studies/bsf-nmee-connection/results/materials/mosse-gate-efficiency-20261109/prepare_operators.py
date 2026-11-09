"""Existing zero-field KS references: linear model gate and planar interface controls.
No SCF/DFT is launched. F_model is an imposed one-electron potential slope, not
a screened external device voltage. All Japan-side processing requires PBS.
"""
from pathlib import Path
import os, socket, sys, json
import numpy as np
from ase.io import read

if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root = Path.cwd()
source = Path('/work/sunq/xjq/mosse-strain-20261108/e000')
parent = Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0] = [str(parent), str(parent/'pawlib')]
from paw import pawpotcar, nonlq
from vaspwfc import vaspwfc

def plane(path):
    with path.open() as f:
        header = [next(f) for _ in range(8)]
        for _ in range(sum(map(int, header[6].split()))): next(f)
        nx, ny, nz = map(int, next(s for s in f if s.strip()).split())
        vals = []
        while len(vals) < nx*ny*nz: vals.extend(map(float, next(f).split()))
    assert len(vals) == nx*ny*nz
    return np.array(vals).reshape(nz, ny, nx).mean(axis=(1, 2))

def shift(v, fraction):
    return np.fft.ifft(np.fft.fft(v)*np.exp(2j*np.pi*np.fft.fftfreq(len(v))*len(v)*fraction)).real

def smooth(v, length, sigma=.20):
    return np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(len(v))*len(v)/length*sigma)**2)).real

def matrix_profiles(seed, profiles):
    atoms = read(seed/'POSCAR')
    wave = vaspwfc(str(seed/'WAVECAR'), lsorbit=True)
    pots = [pawpotcar(s) for s in (seed/'POTCAR').read_text().split('End of Dataset')[:-1]]
    nz = profiles.shape[1]; frac = np.arange(nz)/nz
    fft = np.fft.ifft(profiles, axis=1); nb = wave._nbands; matrices = []
    for ik, k in enumerate(wave._kvecs):
        g = wave.gvectors(ik+1)
        w = np.stack([wave.readBandCoeff(ikpt=ik+1, iband=b+1).reshape(2,-1) for b in range(nb)], axis=-1)
        pr = nonlq(atoms, wave._encut, potcar=pots, k=k)
        out = np.zeros((len(profiles),nb,nb), complex)
        _, groups = np.unique(g[:,:2], axis=0, return_inverse=True)
        for group in np.unique(groups):
            ids = np.flatnonzero(groups == group); gz = g[ids,2]
            block = fft[:,(gz[None,:]-gz[:,None]) % nz]
            for spin in (0,1): out += w[spin,ids].conj().T @ (block @ w[spin,ids])
        for ia, typ in enumerate(pr.element_idx):
            op = pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia]
            betas = [op@ww for ww in w]; q = pots[typ].get_Qij()
            values = np.array([np.interp(atoms.get_scaled_positions()[ia,2], frac, v) for v in profiles])
            out += values[:,None,None]*sum(b.conj().T@q@b for b in betas)
        matrices.append((out+out.conj().transpose(0,2,1))/2)
        print('OPERATOR', seed.parent.name, ik+1, len(wave._kvecs), flush=True)
    return np.array(matrices).transpose(1,0,2,3)

ab = read(source/'A/seed/POSCAR'); mono = read(source/'mono/seed/POSCAR')
vb = plane(source/'A/seed/LOCPOT'); vm = plane(source/'mono/seed/LOCPOT')
assert len(vb) == len(vm)
L = mono.cell[2,2]; z = (np.arange(len(vm))/len(vm)-.5)*L
mo_z = np.sort(ab.positions[np.array(ab.get_chemical_symbols()) == 'Mo',2])
lower = ab[ab.positions[:,2] < L/2]; upper = ab[ab.positions[:,2] > L/2]
zl = lower.positions[:,2]-mo_z[0]; zu = -(upper.positions[:,2]-mo_z[1])
assert lower.get_chemical_symbols() == upper.get_chemical_symbols() == mono.get_chemical_symbols()
upper_mirror_error = float(abs(zl-zu).max())
assert upper_mirror_error < .001, 'Need an independently matched upper-layer reference'
assert abs((mono.positions[:,2]-mono.positions[0,2])-zl).max() < 1e-6
profiles = []
for center, reflected in [(mo_z[0]/L,False), (mo_z[1]/L,True)]:
    mapped = shift(vb, center-.5)
    if reflected: mapped = mapped[(-np.arange(len(vm))) % len(vm)]
    diff = mapped-vm; diff -= diff[abs(z)<1.8].mean(); diff = smooth(diff,L)
    for end in (3.2,3.8):
        u = np.clip((abs(z)-2.2)/(end-2.2),0,1)
        profiles.append(diff*.5*(1+np.cos(np.pi*u)))
gate = z/10  # electron potential energy +e F_model z; eV per V/nm
op_a = matrix_profiles(source/'A/seed',np.array([gate,np.ones_like(gate)]))
op_m = matrix_profiles(source/'mono/seed',np.array([gate,*profiles,np.ones_like(gate)]))
errors = {'A':float(abs(op_a[-1]-np.eye(op_a.shape[-1])).max()),
          'mono':float(abs(op_m[-1]-np.eye(op_m.shape[-1])).max())}
assert max(errors.values()) < 1e-4
np.savez_compressed(root/'operators.npz', A_gate=op_a[0], mono_gate=op_m[0],
                    lower_reference=op_m[1:3],upper_reference=op_m[3:5],
                    profiles=profiles,z_A=z,Mo_separation_A=mo_z[1]-mo_z[0])
record={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'source':str(source),
        'gate':'linear electron potential +F_model*z/10 eV; imposed model slope in V/nm, no new self-consistent screening',
        'matched_lower':'exact intralayer geometry from relaxed e000 bilayer',
        'upper_reference':'horizontal mirror of lower single-layer reference; actual upper intralayer geometry differs slightly',
        'upper_mirror_geometry_error_A':upper_mirror_error,'constant_operator_error':errors,
        'Mo_separation_A':float(mo_z[1]-mo_z[0]),
        'reference':'planar ionic+Hartree difference, smoothed 0.20 A, cutoffs 3.2/3.8 A; nonunique electrostatic reference',
        'no_new_DFT':True}
(root/'operators.json').write_text(json.dumps(record,indent=2));print('OPERATORS_DONE',record,flush=True)
