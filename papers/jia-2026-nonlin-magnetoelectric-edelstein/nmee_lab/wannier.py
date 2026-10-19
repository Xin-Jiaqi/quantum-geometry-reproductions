"""Strict 2D Wannier90 HR adapter for the vetted nondegenerate clean kernel.

H(k)=sum_R exp(i*k.R_nm) H_R/deg_R in a *declared* periodic-Wannier gauge.
The orbital-position connection is set to zero ONLY with explicit opt-in.
For quantitative real-material spin-electric response supply the covariant position
connection and consistently projected spin/Zeeman operators; this module does
not reconstruct them from hr.dat and must not be presented as a DFT-ready solver.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
import numpy as np
from .clean_core.core import BlochModel, layer_spin


@dataclass(frozen=True)
class WannierTB2D:
    translations: np.ndarray  # nr x 2, integer lattice coordinates
    hopping: np.ndarray       # nr x nw x nw, already divided by degeneracy
    lattice_nm: np.ndarray    # 2 x 2; real-space lattice vectors as ROWS
    label: str = 'wannier-model'

    def __post_init__(self):
        r=np.asarray(self.translations)
        h=np.asarray(self.hopping,complex)
        a=np.asarray(self.lattice_nm,float)
        if r.ndim!=2 or r.shape[1]!=2 or not np.allclose(r,np.rint(r)):
            raise ValueError('translations must have integer shape (NR,2)')
        if h.ndim!=3 or h.shape[0]!=len(r) or h.shape[1]!=h.shape[2] or not len(r):
            raise ValueError('hopping shape must be (NR,NW,NW)')
        if a.shape!=(2,2) or abs(np.linalg.det(a))<1e-10:
            raise ValueError('lattice_nm must be a nonsingular 2x2 matrix')
        if not (np.isfinite(h).all() and np.isfinite(a).all()):
            raise ValueError('non-finite hopping/lattice')
        lookup={tuple(x.astype(int)):j for j,x in enumerate(r)}
        if len(lookup)!=len(r):raise ValueError('duplicate R translations')
        for x,j in lookup.items():
            partner=lookup.get(tuple(-np.asarray(x)))
            if partner is None:
                if np.max(abs(h[j]))>1e-12:raise ValueError(f'Hermitian R partner missing for {x}')
            elif not np.allclose(h[j],h[partner].conj().T,atol=2e-9,rtol=2e-9):
                raise ValueError(f'H(R) != H(-R)^dagger at {x}')
        object.__setattr__(self,'translations',r.astype(int))
        object.__setattr__(self,'hopping',h)
        object.__setattr__(self,'lattice_nm',a)

    @property
    def nwann(self):return self.hopping.shape[1]

    @property
    def reciprocal_nm(self):return 2*np.pi*np.linalg.inv(self.lattice_nm).T

    def matrices(self,k_nm,beta=None):
        k=np.asarray(k_nm,float)
        if k.ndim!=2 or k.shape[1]!=2:raise ValueError('k must be (Nk,2) in nm^-1')
        real_R_nm=self.translations@self.lattice_nm
        harmonics=np.exp(1j*k@real_R_nm.T)
        if beta is not None:
            if beta not in (0,1):raise ValueError('beta must be 0/1')
            harmonics*=1j*real_R_nm[:,beta][None,:]
        mat=np.einsum('kr,rij->kij',harmonics,self.hopping,optimize=True)
        return mat


def read_wannier90_hr(path, lattice_nm, *,label='wannier90_hr.dat'):
    """Read conventional Wannier90 hr.dat with 3D R; permit only R_z=0 slabs.

    Original hr.dat rows give H(R) in eV and a per-R degeneracy list.
    lattice_nm MUST be provided separately: Wannier90 hr.dat does not contain it.
    2D systems with physical R_z hopping require a different, 3D adapter.
    """
    lines=Path(path).read_text().splitlines()
    if len(lines)<4:raise ValueError('truncated hr.dat')
    n=int(lines[1].strip()); nr=int(lines[2].strip())
    if n<=0 or nr<=0:raise ValueError('nonpositive dimensions')
    values=[]; pos=3
    while len(values)<nr:
        if pos>=len(lines):raise ValueError('truncated degeneracy table')
        values.extend(int(x) for x in lines[pos].split());pos+=1
    if len(values)!=nr or min(values)<=0:raise ValueError('bad per-R degeneracies')
    content=[line for line in lines[pos:] if line.strip()]
    if len(content)!=nr*n*n:raise ValueError('missing/extra matrix elements')
    Rs=[];Hs=[]
    for ir in range(nr):
        mat=np.zeros((n,n),complex); visited=set(); block=content[ir*n*n:(ir+1)*n*n]
        Rblock=None
        for line in block:
            tok=line.split()
            if len(tok)!=7:raise ValueError('hr data row must have seven columns')
            rx,ry,rz,i,j=map(int,tok[:5])
            if rz!=0:raise ValueError('R_z != 0: not a strictly 2D Wannier HR')
            R=(rx,ry)
            if Rblock is None:Rblock=R
            if R!=Rblock:raise ValueError('unexpected R order within block')
            if not (1<=i<=n and 1<=j<=n) or (i,j) in visited:
                raise ValueError('orbital index invalid or duplicate')
            visited.add((i,j));mat[i-1,j-1]=float(tok[5])+1j*float(tok[6])
        Rs.append(Rblock);Hs.append(mat/values[ir])
    return WannierTB2D(np.array(Rs),np.array(Hs),np.asarray(lattice_nm,float),label=label)


def projector_from_indices(nwann, indices):
    values=[int(i) for i in indices]
    if len(values)!=len(set(values)) or any(i<0 or i>=nwann for i in values):
        raise ValueError('layer indices must be distinct, zero-based, within basis')
    p=np.zeros((nwann,nwann),complex);p[values,values]=1.
    return p


def to_clean_model(tb:WannierTB2D, *, spin_onsite:Mapping[str,np.ndarray],
                   layer_indices:Mapping[str,list[int]]|None=None,g_factor=10.,
                   acknowledge_zero_connection=False):
    """Opt-in exploratory mapping into the established nondegenerate kernel.

    A real Wannier calculation with non-trivial orbital Berry connection or
    non-local spin/position matrices must NOT use this approximation as final
    quantitative material physics. Degeneracy remains rejected by the core.
    """
    if not acknowledge_zero_connection:
        raise ValueError('Explicitly acknowledge the zero-connection / on-site-spin approximation; '
                         'Wannier hr.dat alone is insufficient for physical NMEE')
    n=tb.nwann
    mats={}
    for key in ('x','y','z'):
        if key not in spin_onsite:raise ValueError('need all three spin Pauli matrices, same Wannier basis')
        s=np.asarray(spin_onsite[key],complex)
        if s.shape!=(n,n) or not np.allclose(s,s.conj().T,atol=1e-10):
            raise ValueError('invalid Hermitian spin matrix')
        mats[key]=s
    if layer_indices:
        allidx=[]
        for layer,idx in layer_indices.items():
            p=projector_from_indices(n,idx)
            allidx+=list(idx)
            for component in ('x','y','z'):
                mats[f'{layer}_{component}']=layer_spin(mats[component],p)
        if sorted(allidx)!=list(range(n)):
            raise ValueError('layer indices must partition the entire Wannier basis exactly')
    callbacks={label:(lambda k,m=m:np.broadcast_to(m,(len(k),n,n)).copy()) for label,m in mats.items()}
    return BlochModel(hamiltonian=lambda k:tb.matrices(k),
                      velocity=lambda k,b:tb.matrices(k,b),spin=callbacks,
                      reciprocal=tb.reciprocal_nm,g_factor=g_factor,name=tb.label)
