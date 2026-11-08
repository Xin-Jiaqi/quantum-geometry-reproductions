"""Lower-layer smooth screened electrostatic reference, matched strained geometry.
Not an additive or unique separation of all interface effects. PBS only.
"""
from pathlib import Path
import json,os,socket,sys
import numpy as np
from ase.io import read
from analyze import Response
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute only')
root=Path.cwd();parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104');sys.path[:0]=[str(parent),str(parent/'pawlib')]
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
def shift(v,frac):return np.fft.ifft(np.fft.fft(v)*np.exp(2j*np.pi*np.fft.fftfreq(len(v))*len(v)*frac)).real
def smooth(v,sigma,L):return np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(len(v))*len(v)/L*sigma)**2)).real
base=root/'mono/seed';ab=read(root/'A/seed/POSCAR');at=read(base/'POSCAR');assert np.max(abs(ab.cell.array-at.cell.array))<1e-8
vm=plane(base/'LOCPOT');vb=plane(root/'A/seed/LOCPOT');assert len(vm)==len(vb)
L=at.cell[2,2];nz=len(vm);frac=np.arange(nz)/nz;z=(frac-.5)*L
center=min(ab.positions[[x.symbol=='Mo' for x in ab],2])/L
# The isolated layer is centered by its Mo plane and has identical intralayer geometry.
lo=ab[ab.positions[:,2]<L/2];lo.positions[:,2]+=(.5-center)*L
assert np.max(abs(lo.positions-at.positions))<1e-6
profiles=[];diff=shift(vb,center-.5)-vm;diff-=diff[abs(z)<1.8].mean();diff=smooth(diff,.20,L)
for end in (3.2,3.8):
 u=np.clip((abs(z)-2.2)/(end-2.2),0,1);profiles.append(diff*.5*(1+np.cos(np.pi*u)))
profiles=np.array(profiles);vfft=np.fft.ifft(profiles,axis=1);wave=vaspwfc(str(base/'WAVECAR'),lsorbit=True);nb=wave._nbands
pots=[pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]];operators=[]
for ik,k in enumerate(wave._kvecs):
 g=wave.gvectors(ik+1);w=np.stack([wave.readBandCoeff(ikpt=ik+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1);pr=nonlq(at,wave._encut,potcar=pots,k=k)
 out=np.zeros((2,nb,nb),complex);_,groups=np.unique(g[:,:2],axis=0,return_inverse=True)
 for group in np.unique(groups):
  ids=np.where(groups==group)[0];gz=g[ids,2];block=vfft[:,(gz[None,:]-gz[:,None])%nz]
  for spin in (0,1):out+=w[spin,ids].conj().T@(block@w[spin,ids])
 for ia,typ in enumerate(pr.element_idx):
  op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia];betas=[op@ww for ww in w];q=pots[typ].get_Qij();values=np.array([np.interp(at.get_scaled_positions()[ia,2],frac,v) for v in profiles])
  out+=values[:,None,None]*sum(b.conj().T@q@b for b in betas)
 operators.append((out+out.conj().transpose(0,2,1))/2)
 print('REFERENCE_OPERATOR',ik+1,len(wave._kvecs),flush=True)
operators=np.array(operators).transpose(1,0,2,3);np.savez_compressed(root/'reference_operators.npz',operators=operators,profiles=profiles,z_A=z)
d=np.load(root/'mono/cache.npz');result={'job':os.environ['PBS_JOBID'],'node':socket.gethostname(),'matched_intralayer_geometry':True,'units':'1e-9 (nm V T)^-1','scope':__doc__,'rows':[]}
for nb in (48,56,64):
 r=Response(d,nb)
 for i,end in enumerate((3.2,3.8)):
  ans=r.response(r.h+operators[i,:,:nb,:nb]);result['rows'].append(dict(nb=nb,taper_end_A=end,C_lower=ans['C_total_lower_upper'][0],**ans))
(root/'reference.json').write_text(json.dumps(result,indent=2));np.savez_compressed(root/'reference_profiles.npz',profiles=profiles,z_A=z)
print('REFERENCE_DONE',flush=True)
