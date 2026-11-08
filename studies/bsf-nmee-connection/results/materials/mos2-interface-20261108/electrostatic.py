"""Test smooth interface electrostatics on isolated MoS2, without cross-layer states.

Planar ionic+Hartree potential only; no XC/nonlocal potential, local-field Fourier
components, or self-consistent screening. PAW augmentation uses V at each atom.
The interface potential is smoothly cut outside the isolated layer, avoiding the
other layer's ionic well. Two cutoff radii test this external-embedding choice.
"""
from pathlib import Path
import os,socket,sys,json
import numpy as np
from ase.io import read
from analyze import Response
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
root=Path.cwd();parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0]=[str(parent),str(parent/'pawlib')]
from vaspwfc import vaspwfc
from paw import pawpotcar,nonlq

def plane(path):
    with path.open() as f:
        head=[next(f) for _ in range(8)]
        for _ in range(sum(map(int,head[6].split()))):next(f)
        dims=list(map(int,next(s for s in f if s.strip()).split()));nx,ny,nz=dims
        vals=[]
        while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
    assert len(vals)==nx*ny*nz
    return np.array(vals).reshape(nz,ny,nx).mean(axis=(1,2))
def shift(v,frac):
    return np.fft.ifft(np.fft.fft(v)*np.exp(2j*np.pi*np.fft.fftfreq(len(v))*len(v)*frac)).real
def smooth(v,sigma,L):
    return np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(len(v))*len(v)/L*sigma)**2)).real
mono=root/'mono'/'seed';bi=root/'A'/'seed';at=read(mono/'POSCAR');ab=read(bi/'POSCAR')
vm=plane(mono/'LOCPOT');vb=plane(bi/'LOCPOT');assert len(vm)==len(vb)
L=at.cell[2,2];nz=len(vm);frac=np.arange(nz)/nz;zc=(frac-.5)*L
mo=ab.positions[[i for i,x in enumerate(ab) if x.symbol=='Mo'],2];centers=sorted(mo/L)
# Difference of ionic+Hartree potentials, shifted to the isolated layer center.
profiles={};meta={}
for layer,center in zip(('lower','upper'),centers):
    diff=shift(vb,center-.5)-vm
    # Global energy offset is arbitrary; remove its mean within the constituent.
    diff-=diff[abs(zc)<1.8].mean();diff=smooth(diff,.20,L)
    for end in (3.2,3.8):
        start=2.2;window=np.ones(nz);u=np.clip((abs(zc)-start)/(end-start),0,1);window=.5*(1+np.cos(np.pi*u))
        v=diff*window;profiles[f'{layer}_{end}_full']=v
        profiles[f'{layer}_{end}_odd']=(v-v[(-np.arange(nz))%nz])/2
    sel=abs(zc)<1.8;slope=float(np.polyfit(zc[sel],diff[sel],1)[0])
    meta[layer]={'linear_fit_physical_Vnm':10*slope,'deltaV_S_planes_eV':float(np.interp(1.56,zc,diff)-np.interp(-1.56,zc,diff)),
                 'scope':'plane-average interface potential, independently calculated, not fitted to response'}
# Independent calibration against the prior isolated-layer gate calculation.
for end in (3.2,3.8):
    u=np.clip((abs(zc)-2.2)/(end-2.2),0,1)
    profiles[f'uniform_F+0.25_{end}']=.025*zc*.5*(1+np.cos(np.pi*u))
# Purely interface-induced Hartree redistribution; ions cancel in continuous space.
redis=vb-sum(shift(vm,.5-c) for c in centers);redis-=redis.mean()
np.savez_compressed(root/'potential_profiles.npz',z_A=zc,mono=vm,bilayer=vb,redistribution=smooth(redis,.20,L),**profiles)
(root/'electrostatic_profiles.json').write_text(json.dumps(meta,indent=2))
wave=vaspwfc(str(mono/'WAVECAR'),lsorbit=True);cache=np.load(root/'mono'/'cache.npz');r=Response(cache)
pots=[pawpotcar(s) for s in (mono/'POTCAR').read_text().split('End of Dataset')[:-1]]
nb=wave._nbands; names=list(profiles);vs=np.array([profiles[s] for s in names]);vfft=np.fft.ifft(vs,axis=1)
operators=[]
for ik,kv in enumerate(wave._kvecs):
    g=wave.gvectors(ik+1);w=np.stack([wave.readBandCoeff(ikpt=ik+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1)
    pr=nonlq(at,wave._encut,potcar=pots,k=kv);out=np.zeros((len(vs),nb,nb),complex)
    _,groups=np.unique(g[:,:2],axis=0,return_inverse=True)
    for group in np.unique(groups):
        ids=np.where(groups==group)[0];gz=g[ids,2];block=vfft[:,(gz[None,:]-gz[:,None])%nz]
        for spin in (0,1):out+=w[spin,ids].conj().T@(block@w[spin,ids])
    for ia,typ in enumerate(pr.element_idx):
        op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia];betas=[op@v for v in w];q=pots[typ].get_Qij()
        value=np.array([np.interp(at.get_scaled_positions()[ia,2],frac,v) for v in vs])
        out+=value[:,None,None]*sum(b.conj().T@q@b for b in betas)
    operators.append((out+out.conj().transpose(0,2,1))/2)
    print('ELECTROSTATIC_EXPORT',ik+1,len(wave._kvecs),flush=True)
operators=np.array(operators).transpose(1,0,2,3)
np.savez_compressed(root/'electrostatic_operators.npz',names=names,operators=operators)
result={'scope':__doc__,'profiles':meta,'units':'1e-9 (nm V T)^-1','baseline':r.response(),'controls':{}}
for name,v in zip(names,operators):
    result['controls'][name]=r.response(r.h+v)
    (root/'electrostatic.json').write_text(json.dumps(result,indent=2))
print('ELECTROSTATIC_DONE',result,flush=True)
