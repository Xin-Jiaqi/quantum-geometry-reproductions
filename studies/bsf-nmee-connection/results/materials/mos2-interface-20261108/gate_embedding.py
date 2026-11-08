"""Calibrate planar electrostatic embedding using one already-tested gate value."""
from pathlib import Path
import os,socket,sys,json,numpy as np
from ase.io import read
from analyze import Response
if not os.getenv('PBS_JOBID') or socket.gethostname().startswith(('super','cms-ssh')):raise SystemExit('PBS compute required')
parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104');sys.path[:0]=[str(parent),str(parent/'pawlib')]
from vaspwfc import vaspwfc
from paw import pawpotcar,nonlq
root=Path.cwd();base=root/'mono/seed';gated=root/'mono_F+0.250/seed';at=read(base/'POSCAR')
def plane(path):
 with path.open() as f:
  h=[next(f) for _ in range(8)]
  for _ in range(sum(map(int,h[6].split()))):next(f)
  nx,ny,nz=map(int,next(s for s in f if s.strip()).split());vals=[]
  while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
 return np.array(vals).reshape(nz,ny,nx).mean(axis=(1,2))
v=plane(gated/'LOCPOT')-plane(base/'LOCPOT');nz=len(v);frac=np.arange(nz)/nz;L=at.cell[2,2];z=(frac-.5)*L
v-=v[abs(z)<1.8].mean();v=np.fft.ifft(np.fft.fft(v)*np.exp(-.5*(2*np.pi*np.fft.fftfreq(nz)*nz/L*.20)**2)).real
vs=[v];names=['whole_cell']
for end in (3.2,3.8):
 window=.5*(1+np.cos(np.pi*np.clip((abs(z)-2.2)/(end-2.2),0,1)))
 vs.append(v*window);names.append(f'cutoff_{end}')
vs=np.array(vs);fft=np.fft.ifft(vs,axis=1);wave=vaspwfc(str(base/'WAVECAR'),lsorbit=True);nb=wave._nbands
pots=[pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]];operators=[]
for ik,k in enumerate(wave._kvecs):
 g=wave.gvectors(ik+1);w=np.stack([wave.readBandCoeff(ikpt=ik+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1)
 pr=nonlq(at,wave._encut,potcar=pots,k=k);out=np.zeros((len(vs),nb,nb),complex);_,groups=np.unique(g[:,:2],axis=0,return_inverse=True)
 for group in np.unique(groups):
  ids=np.where(groups==group)[0];gz=g[ids,2];block=fft[:,(gz[None,:]-gz[:,None])%nz]
  for s in (0,1):out+=w[s,ids].conj().T@(block@w[s,ids])
 for ia,typ in enumerate(pr.element_idx):
  op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia];betas=[op@ww for ww in w];q=pots[typ].get_Qij()
  values=np.array([np.interp(at.get_scaled_positions()[ia,2],frac,vv) for vv in vs]);out+=values[:,None,None]*sum(b.conj().T@q@b for b in betas)
 operators.append((out+out.conj().transpose(0,2,1))/2)
operators=np.array(operators).transpose(1,0,2,3);d=np.load(root/'mono/cache.npz')
result={'scope':'Independent self-consistent Fz=+0.25 V/nm potential on zero-field monolayer; planar electrostatics only',
 'physical_external_Vnm':.25,'vacuum_slope_Vnm':float(10*np.polyfit(z[(z>5)&(z<9)],v[(z>5)&(z<9)],1)[0]),
 'layer_linear_fit_Vnm':float(10*np.polyfit(z[abs(z)<1.8],v[abs(z)<1.8],1)[0]),'cutoffs':{}}
for nb in (40,48,56,64):
 r=Response(d,nb);result['cutoffs'][str(nb)]={name:r.response(r.h+op[:,:nb,:nb]) for name,op in zip(names,operators)}
(root/'gate_embedding.json').write_text(json.dumps(result,indent=2));np.savez_compressed(root/'gate_potential.npz',z_A=z,profiles=vs,names=names)
print('GATE_EMBEDDING_DONE',result,flush=True)
