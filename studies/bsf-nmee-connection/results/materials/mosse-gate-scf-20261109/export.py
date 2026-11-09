"""Full-KS links, spin/layer operators and atomic trials; no downfolding of response."""
from pathlib import Path
import os,socket,sys,json
import numpy as np
from ase.io import read
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute allocation required')
parent=Path('/work/sunq/xjq/mosse-7269-mechanism-20261104')
sys.path[:0]=[str(parent),str(parent/'pawlib')]
from paw import pawpotcar,nonlq
from vaspwfc import vaspwfc
from spin_transitions import spin_and_metric,fft_z
root=Path.cwd(); label=os.environ['STRUCTURE'];d=root/label/'seed'
a=read(d/'POSCAR'); wave=vaspwfc(str(d/'WAVECAR'),lsorbit=True)
pots=[pawpotcar(s) for s in (d/'POTCAR').read_text().split('End of Dataset')[:-1]]
nb=wave._nbands;k=wave._kvecs;nk=len(k);n=round(nk**.5);assert n*n==nk
symbols=list(dict.fromkeys(a.get_chemical_symbols())); ncore=10*sum(x.symbol=='Mo' for x in a)
nocc=int((d/'EIGENVAL').read_text().splitlines()[5].split()[0]); assert nocc==24*(2 if label=='A' else 1)
z=a.get_scaled_positions()[:,2];nz=fft_z(d/'CHGCAR');mask=np.zeros(nz);mask[1:nz//2]=1;mask[[0,nz//2]]=.5
atom_bottom=z<.5
orb_bottom=np.repeat(atom_bottom,[10 if x.symbol=='Mo' else 6 for x in a])
terms=[]
for atom in a:
    l=2 if atom.symbol=='Mo' else 1;pp=pots[symbols.index(atom.symbol)]
    first=next(q for q,ll,m in pp.ilm if ll==l)
    def radial(q):return pp.radial_simp_int(pp.paw_ae_wfc[first]*pp.paw_ae_wfc[q],inside_rcomp=True)
    norm=np.sqrt(radial(first))
    terms.append(np.array([[radial(q)/norm if ll==l and mm==m else 0. for q,ll,mm in pp.ilm] for m in range(-l,l+1)]))
steps=[(1,0),(-1,0),(0,1),(0,-1)]
links=np.zeros((nk,4,nb,nb),complex);neighbors=np.full((nk,4),-1)
with (d/'wannier90.mmn').open() as f:
    next(f);bands,nn,nh=map(int,next(f).split());assert (bands,nn)==(nb,nk)
    for _ in range(nk*nh):
        i,j,*g=map(int,next(f).split());i-=1;j-=1
        block=np.fromstring(' '.join(next(f) for _ in range(nb*nb)),sep=' ').reshape(-1,2)
        vec=k[j]+g-k[i];b=tuple(np.rint(vec[:2]*n).astype(int))
        if abs(vec[2])<1e-5 and b in steps:
            q=steps.index(b);neighbors[i,q]=j;links[i,q]=(block[:,0]+1j*block[:,1]).reshape(nb,nb).T
assert np.all(neighbors>=0)
# A link and its reverse must be adjoints, including reciprocal-lattice boundary phases.
reverse_error=max(np.max(abs(links[i,t]-links[j,t^1].conj().T)) for i in range(nk) for t,j in enumerate(neighbors[i]))
if reverse_error>1e-6:raise RuntimeError(f'nonreciprocal link {reverse_error}')
# Independent check of the eigenvector gauge shared by WAVECAR and MMN.
i=next(q for q in range(nk) if abs(k[q]).max()<1e-8);j=neighbors[i,0]
def states(ix):
    return wave.gvectors(ix+1),np.stack([wave.readBandCoeff(ikpt=ix+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1)
gi,wi=states(i);gj,wj=states(j);sh=np.rint(k[i]+[1/n,0,0]-k[j]).astype(int)
jmap={tuple(g):ii for ii,g in enumerate(gj)};ii=[q for q,g in enumerate(gi) if tuple(g+sh) in jmap];jj=[jmap[tuple(gi[q]+sh)] for q in ii]
approx=sum(wi[s,ii].conj().T@wj[s,jj] for s in (0,1));pr=[nonlq(a,wave._encut,potcar=pots,k=k[q]) for q in (i,j)]
for ia,typ in enumerate(pr[0].element_idx):
    pp=[v.qproj[typ]*v.cqfak[typ][:,None]*v.crexp[None,:,ia] for v in pr]
    phase=np.exp(-2j*np.pi*sh@a.get_scaled_positions()[ia])
    for spin in (0,1):approx+=phase*(pp[0]@wi[spin]).conj().T@pots[typ].get_Qij()@(pp[1]@wj[spin])
sl=slice(ncore,nocc+8);gauge=float(np.linalg.norm((approx-links[i,0])[sl,sl])/np.linalg.norm(links[i,0][sl,sl]))
if gauge>.08:raise RuntimeError(f'Wave/link gauge mismatch {gauge}')
coeff=[];S=[];ST=[];Z=[];P=[];metric=[]
for i,kv in enumerate(k):
    g=wave.gvectors(i+1);w=np.stack([wave.readBandCoeff(ikpt=i+1,iband=b+1).reshape(2,-1) for b in range(nb)],axis=-1)
    pr=nonlq(a,wave._encut,potcar=pots,k=kv);assert np.array_equal(g,pr.Gvec)
    ov,s=spin_and_metric(w,pr,g,mask,atom_bottom);err=float(abs(ov-np.eye(nb)).max());metric.append(err)
    if err>1e-4:raise RuntimeError(f'bad PAW norm {err}')
    sx=w[0].conj().T@w[1]+w[1].conj().T@w[0];sz=w[0].conj().T@w[0]-w[1].conj().T@w[1]
    cb=[]
    for ia,typ in enumerate(pr.element_idx):
        op=pr.qproj[typ]*pr.cqfak[typ][:,None]*pr.crexp[None,:,ia]
        b0,b1=[op@v for v in w];q=pots[typ].get_Qij()
        sx+=b0.conj().T@q@b1+b1.conj().T@q@b0;sz+=b0.conj().T@q@b0-b1.conj().T@q@b1
        phase=np.exp(2j*np.pi*(kv@a.get_scaled_positions()[ia]))
        cb.append(np.stack([phase*terms[ia]@b[:,ncore:] for b in (b0,b1)],axis=1).reshape(-1,nb-ncore))
    combos=np.zeros_like(w);combos[0]=combos[1]=w[0]/np.sqrt(2)
    _,p0=spin_and_metric(combos,pr,g,mask,atom_bottom)
    combos[0]=combos[1]=w[1]/np.sqrt(2);_,p1=spin_and_metric(combos,pr,g,mask,atom_bottom)
    S.append(s);ST.append(sx);Z.append(sz);P.append(p0+p1);coeff.append(np.concatenate(cb))
    print('EXPORT',label,i+1,nk,flush=True)
energy=wave._bands[0];dest=root/label/'cache.npz'
np.savez_compressed(dest,energies=energy,Z=Z,spin=S,spin_total=ST,P=P,links=links,neighbors=neighbors,
    k=k,cell=a.cell.array,pos=a.get_scaled_positions(),coeff=coeff,nocc=nocc,ncore=ncore,orb_bottom=orb_bottom)
record={'label':label,'bands':nb,'nocc':nocc,'nk':nk,'metric_max':max(metric),'link_reverse_error':float(reverse_error),'wave_link_gauge_error':gauge,
        'gap_eV':float(energy[:,nocc].min()-energy[:,nocc-1].max()),'job':os.environ['PBS_JOBID'],'node':socket.gethostname()}
(root/label/'export.json').write_text(json.dumps(record,indent=2));print('EXPORTED',record,flush=True)
