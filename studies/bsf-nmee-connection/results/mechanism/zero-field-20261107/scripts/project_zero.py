"""Localized d+p trial overlaps for the 7269 model. Run only under PBS."""
import base64
import os
from pathlib import Path
import socket
import sys
import numpy as np
from ase.io import read

if not os.environ.get('PBS_JOBID') or not os.environ.get('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('compute allocation required')
root=Path.cwd()
sys.path[:0]=[str(root.parent),str(root.parent/'pawlib')]
from paw import pawpotcar,nonlq
from vaspwfc import vaspwfc
tag=os.environ['MAT']
base=root/'seed'
atoms=read(base/'POSCAR')
wave=vaspwfc(str(base/'WAVECAR'),lsorbit=True)
pots=[pawpotcar(s) for s in (base/'POTCAR').read_text().split('End of Dataset')[:-1]]
# Use one normalized AE radial trial per l, contract *all* projector channels for that l.
# Coefficients from nonlq carry e^iG.tau. Restore e^ik.tau for periodic orbital gauge.
terms=[]
for ia,atom in enumerate(atoms):
    l=2 if atom.symbol=='Mo' else 1
    pp=pots[['Mo','S','Se'].index(atom.symbol)]
    first=next(n for n,ll,m in pp.ilm if ll==l)
    def integral(n):
        return pp.radial_simp_int(pp.paw_ae_wfc[first]*pp.paw_ae_wfc[n],inside_rcomp=True)
    norm=np.sqrt(integral(first))
    terms.append(np.array([[integral(n)/norm if ll==l and mm==m else 0.
                            for n,ll,mm in pp.ilm] for m in range(-l,l+1)]))
coeff=[]
for ik,k in enumerate(wave._kvecs):
    projector=nonlq(atoms,wave._encut,potcar=pots,k=k)
    assert np.array_equal(projector.Gvec,wave.gvectors(ik+1))
    states=np.stack([wave.readBandCoeff(ikpt=ik+1,iband=b+1).reshape(2,-1) for b in range(20,112)],axis=-1)
    blocks=[]
    for ia,typ in enumerate(projector.element_idx):
        op=terms[ia]@(projector.qproj[typ]*projector.cqfak[typ][:,None])
        op*=projector.crexp[None,:,ia]
        phase=np.exp(2j*np.pi*(k@atoms.get_scaled_positions()[ia]))
        blocks.append(np.stack([phase*(op@spin) for spin in states],axis=1).reshape(-1,92))
    coeff.append(np.concatenate(blocks))
out=root/f'dp-{tag}.npz'
np.savez_compressed(out,coeff=np.array(coeff,dtype=np.complex64),energy=wave._bands[0,:,20:112],
                    k=wave._kvecs,cell=atoms.cell.array,pos=atoms.get_scaled_positions())
