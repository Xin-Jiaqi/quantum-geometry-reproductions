"""Discrete Berry electric coupling using first-principles overlaps on the original 6x6 mesh."""
from pathlib import Path
import json,os,sys,socket
import numpy as np
from scipy.linalg import null_space
from ase.io import read
if not os.environ.get('PBS_JOBID') or not os.environ.get('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('compute allocation required')
root=Path.cwd()
sys.path[:0]=[str(root.parent),str(root.parent/'pawlib')]
from paw import pawpotcar,nonlq
from vaspwfc import vaspwfc
from model import polar_cols,localized_model,BOTTOM
from spin_transitions import spin_and_metric,fft_z
tag=os.environ['MAT'];src=np.load(f'dp-{tag}.npz');k=src['k'];cell=src['cell']
d=root/'seed';waves=vaspwfc(str(d/'WAVECAR'),lsorbit=True)
a=read(d/'POSCAR');pots=[pawpotcar(s) for s in (d/'POTCAR').read_text().split('End of Dataset')[:-1]]
z=a.get_scaled_positions()[:,2];cut=np.sort(z)[2:4].mean();nz=fft_z(d/'CHGCAR')
mask=((np.arange(nz)/nz-cut)%1>=.5).astype(float)
nk=len(k);mesh=int(round(np.sqrt(nk)));steps=[(1,0),(-1,0),(0,1),(0,-1)]
links=np.empty((nk,4,112,112),complex);neighbors=np.full((nk,4),-1)
with (root/'seed'/'wannier90.mmn').open() as f:
    next(f);nb,nn,nbh=map(int,next(f).split());assert (nb,nn)==(112,nk)
    for _ in range(nk*nbh):
        i,j,*g=map(int,next(f).split());i-=1;j-=1
        block=np.fromstring(' '.join(next(f) for _ in range(nb*nb)),sep=' ').reshape(-1,2)
        vec=k[j]+g-k[i];b=tuple(np.rint(vec[:2]*mesh).astype(int))
        if abs(vec[2])<1e-5 and b in steps:
            axis=steps.index(b);neighbors[i,axis]=j;links[i,axis]=(block[:,0]+1j*block[:,1]).reshape(nb,nb).T
assert np.all(neighbors>=0)
def read_w(i):
    g=waves.gvectors(i+1)
    return g,np.stack([waves.readBandCoeff(ikpt=i+1,iband=j+1).reshape(2,-1) for j in range(112)],axis=-1)
# Independent gauge sanity test: first neighboring link, using pseudo plane waves and
# the leading PAW augmentation form factor. Small nonzero differences are expected.
i=next(ix for ix in range(nk) if abs(k[ix]).max()<1e-5);j=neighbors[i,0]
gi,wi=read_w(i);gj,wj=read_w(j)
shift=np.rint(k[i]+np.array([1/mesh,0,0])-k[j]).astype(int)
jl={tuple(row):idx for idx,row in enumerate(gj)}
ix=[idx for idx,row in enumerate(gi) if tuple(row+shift) in jl];jx=[jl[tuple(gi[t]+shift)] for t in ix]
approx=sum(wi[s,ix].conj().T@wj[s,jx] for s in (0,1))
pr=[nonlq(a,waves._encut,potcar=pots,k=k[ix0]) for ix0 in (i,j)]
for ia,typ in enumerate(pr[0].element_idx):
    ops=[obj.qproj[typ]*obj.cqfak[typ][:,None]*obj.crexp[None,:,ia] for obj in pr]
    phase=np.exp(-2j*np.pi*(shift@a.get_scaled_positions()[ia]))
    for s in (0,1):
        approx+=phase*(ops[0]@wi[s]).conj().T@pots[typ].get_Qij()@(ops[1]@wj[s])
q=slice(20,56)
gauge_error=np.linalg.norm((approx-links[i,0])[q,q])/np.linalg.norm(links[i,0][q,q])
if gauge_error>.08: raise RuntimeError(f'Wavefunction and link gauge mismatch {gauge_error}')
def chi(e,u,spin,mats,nocc,pairs=None):
    result=0.
    for i in range(nk):
        occ=u[i,:,:nocc];empty=u[i,:,nocc:];deriv=0.
        for t,j in enumerate(neighbors[i]):
            shifted=mats[i,t]@u[j,:,:nocc]
            dual=np.linalg.solve((occ.conj().T@shifted).T,shifted.T).T
            deriv+=((1 if t%2==0 else -1)*mesh*cell[t//2,0]/(4*np.pi))*dual
        dip=empty.conj().T@(1j*deriv)
        spinvc=occ.conj().T@spin[i]@empty
        p=2*np.real(spinvc*dip.T)/(e[i,:nocc,None]-e[i,None,nocc:])
        result+=p.sum()
        if pairs is not None:pairs.append(p)
    return result/nk
nempty=int(os.environ.get('FROZEN_EMPTY','8'));nf=28+nempty
spin_all=[];lift=[];p_all=[];z_all=[];sx_all=[]
for i in range(nk):
    g,w=read_w(i);pr=nonlq(a,waves._encut,potcar=pots,k=k[i])
    O,S=spin_and_metric(w,pr,g,mask,z<cut);assert abs(O-np.eye(112)).max()<1e-4
    spin_all.append(S)
    zs=w[0].conj().T@w[0]-w[1].conj().T@w[1]
    sx=w[0].conj().T@w[1]+w[1].conj().T@w[0]
    for ia,typ in enumerate(pr.element_idx):
        proj=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia]
        b0,b1=[proj@v for v in w];q=pr.pawpp[typ].get_Qij()
        zs+=b0.conj().T@q@b0-b1.conj().T@q@b1
        sx+=b0.conj().T@q@b1+b1.conj().T@q@b0
    z_all.append((zs+zs.conj().T)/2)
    sx_all.append((sx+sx.conj().T)/2)

    # <P_layer> is obtained by summing the two diagonal spin blocks
    # with the same PAW step-function integration as Sx.
    combos=np.zeros_like(w)
    combos[0]=combos[1]=w[0]/np.sqrt(2)
    _,pup=spin_and_metric(combos,pr,g,mask,z<cut)
    combos[0]=combos[1]=w[1]/np.sqrt(2)
    _,pdown=spin_and_metric(combos,pr,g,mask,z<cut)
    p_all.append(pup+pdown)
    cf=src['coeff'][i].astype(complex);e=src['energy'][i];f,_=polar_cols(cf[:,:nf]);other=null_space(f.conj().T)
    t,_=polar_cols(((other.conj().T@cf[:,nf:])/(1+((e[nf:]-e[nf])/4)**2)).conj().T)
    lift.append(np.vstack([np.zeros((20,44)),f.conj().T,t@other.conj().T]))
lift=np.array(lift)
hg,_=localized_model(src,nempty);h=np.array([hg[tuple(np.rint(mesh*x[:2]).astype(int)%mesh)] for x in k])
# Perturb only the chosen d+p interlayer block in the 112-state DFT space;
# retain its actual Bloch links, other states and full real-space spin readout.
inter=BOTTOM[:,None]!=BOTTOM[None,:]
from model import NAMES
names=np.repeat(NAMES,2);el=np.repeat(['Mo','Mo','S','S','Se','Se'],[10,10,6,6,6,6])
ip=(el=='S')&(names=='pz')
pzpair=inter & (ip[:,None]&ip[None,:]);withpz=inter & (ip[:,None]|ip[None,:])
deep=np.isin(names,['dxz','dyz','px','py'])
deep=inter & ~withpz & (deep[:,None]|deep[None,:])
blocks={'all':inter,'S pz pair':pzpair,'all S pz':withpz,'other d1 pxy':deep,'remaining':inter & ~withpz & ~deep}
En=waves._bands[0];diag=En[:,:,None]*np.eye(112)
spinKS=np.array(spin_all);pKS=np.array(p_all)

Z=np.array(z_all)
np.savez_compressed('full-controls-cache.npz',diag=diag,Z=Z,spin=spinKS,spin_total=np.array(sx_all),P=pKS,links=links,neighbors=neighbors,k=k,cell=cell,coeff=src['coeff'],energy=src['energy'])
Path('full-controls-export.json').write_text(json.dumps({'mesh':mesh,'nk':nk,'bands':112,'source_field_eV':0.,'gauge_error':gauge_error,'node':socket.gethostname(),'job':os.environ['PBS_JOBID']}))
print('FULL_EXPORT_COMPLETE',flush=True)
