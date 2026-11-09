"""Smooth self-consistent layer environments on isolated MoS2; PBS only.

This is a screened electrostatic representation test, not a unique removal of
all effects caused by hybridization. The parent density is the full bilayer's.
"""
from pathlib import Path
import json,os,socket,sys
import numpy as np
from ase.io import read
from analyze import Response
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute allocation required')
root=Path.cwd();old=Path('/work/sunq/xjq/mos2-interface-20261108')
parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0]=[str(parent),str(parent/'pawlib')]
from paw import pawpotcar,nonlq
from vaspwfc import vaspwfc
def plane(path):
    with path.open() as f:
        h=[next(f) for _ in range(8)]
        for _ in range(sum(map(int,h[6].split()))):next(f)
        nx,ny,nz=map(int,next(s for s in f if s.strip()).split());vals=[]
        while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
    assert len(vals)==nx*ny*nz
    return np.array(vals).reshape(nz,ny,nx).mean(axis=(1,2))
def shift(v,frac):
    return np.fft.ifft(np.fft.fft(v)*np.exp(2j*np.pi*np.fft.fftfreq(len(v))*len(v)*frac)).real
def smooth(v,sigma,L):
    return np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(len(v))*len(v)/L*sigma)**2)).real
cases=json.loads((root/'reference_cases.json').read_text())
base=old/'mono/seed';at=read(base/'POSCAR');vm=plane(base/'LOCPOT')
nz=len(vm);L=at.cell[2,2];frac=np.arange(nz)/nz;z=(frac-.5)*L
profiles=[];records=[];profile_meta=[]
for case in cases:
    seed=Path(case['seed']);assert (seed/'finished.json').exists()
    ab=read(seed/'POSCAR');assert np.max(abs(ab.cell.array-at.cell.array))<1e-8
    vb=plane(seed/'LOCPOT');assert len(vb)==nz
    centers=sorted(ab.positions[[i for i,a in enumerate(ab) if a.symbol=='Mo'],2]/L)
    for layer,center in zip(('lower','upper'),centers):
        diff=shift(vb,center-.5)-vm
        diff-=diff[abs(z)<1.8].mean();diff=smooth(diff,.20,L)
        for end in (3.2,3.8):
            u=np.clip((abs(z)-2.2)/(end-2.2),0,1);profile=diff*.5*(1+np.cos(np.pi*u))
            profiles.append(profile);records.append(dict(case=case['case'],layer=layer,taper_end_A=end))
            profile_meta.append({'case':case['case'],'layer':layer,'taper_end_A':end,
                'internal_linear_fit_Vnm':float(10*np.polyfit(z[abs(z)<1.8],diff[abs(z)<1.8],1)[0])})
profiles=np.array(profiles);vfft=np.fft.ifft(profiles,axis=1)
wave=vaspwfc(str(base/'WAVECAR'),lsorbit=True);nb=wave._nbands
pots=[pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]]
operators=[]
for ik,k in enumerate(wave._kvecs):
    g=wave.gvectors(ik+1);w=np.stack([wave.readBandCoeff(ikpt=ik+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1)
    pr=nonlq(at,wave._encut,potcar=pots,k=k)
    out=np.zeros((len(profiles),nb,nb),complex);_,groups=np.unique(g[:,:2],axis=0,return_inverse=True)
    for group in np.unique(groups):
        ids=np.where(groups==group)[0];gz=g[ids,2];block=vfft[:,(gz[None,:]-gz[:,None])%nz]
        for spin in (0,1):out+=w[spin,ids].conj().T@(block@w[spin,ids])
    for ia,typ in enumerate(pr.element_idx):
        op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia];betas=[op@ww for ww in w];q=pots[typ].get_Qij()
        values=np.array([np.interp(at.get_scaled_positions()[ia,2],frac,v) for v in profiles])
        out+=values[:,None,None]*sum(b.conj().T@q@b for b in betas)
    operators.append((out+out.conj().transpose(0,2,1))/2)
    print('REFERENCE_OPERATOR',ik+1,len(wave._kvecs),flush=True)
operators=np.array(operators).transpose(1,0,2,3)
np.savez_compressed(root/'reference_operators.npz',operators=operators,profiles=profiles,z_A=z)
d=np.load(old/'mono/cache.npz')
result={'scope':__doc__,'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),
        'cases':cases,'profile_metadata':profile_meta,'units':'1e-9 (nm V T)^-1, Pauli spin density','rows':[]}
for cutoff in (40,48,56,64):
    r=Response(d,cutoff)
    for case in cases:
        for taper in (3.2,3.8):
            vals={};gaps=[]
            for layer in ('lower','upper'):
                index=records.index(dict(case=case['case'],layer=layer,taper_end_A=taper))
                ans=r.response(r.h+operators[index,:,:cutoff,:cutoff])
                # Each entire isolated monolayer represents one material layer.
                vals[layer]=ans['C_total_lower_upper'][0];gaps.append(ans['gap_eV'])
            result['rows'].append({'case':case['case'],'delta_d_A':case['delta_d_A'],'Fz_Vnm':case['Fz_Vnm'],
                'cutoff':cutoff,'taper_end_A':taper,'C_lower':vals['lower'],'C_upper':vals['upper'],
                'C_total':vals['lower']+vals['upper'],'gap_min_eV':min(gaps)})
    (root/'electrostatic_reference.json').write_text(json.dumps(result,indent=2))
np.savez_compressed(root/'reference_profiles.npz',profiles=profiles,z_A=z)
print('ELECTROSTATIC_REFERENCE_DONE',len(result['rows']),flush=True)
