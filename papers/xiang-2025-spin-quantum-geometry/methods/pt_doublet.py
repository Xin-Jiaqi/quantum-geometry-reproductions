"""PT-symmetric rank-two subspace check, migrated from offline V16.

Synthetic 4-band gamma model, NOT MnBi2Te4 or a material prediction.
The projector resolvent is independent of the band/orbital Liouville method.
"""
from __future__ import annotations
import numpy as np

Sx=np.array([[0,1],[1,0]],complex)
Sy=np.array([[0,-1j],[1j,0]],complex)
Sz=np.diag([1.,-1.]).astype(complex)
I=np.eye(2,dtype=complex)
I4=np.eye(4,dtype=complex)
G=np.stack([np.kron(Sx,I),np.kron(Sy,I),np.kron(Sz,Sx),
            np.kron(Sz,Sy),np.kron(Sz,Sz)])
PT=np.kron(Sx,1j*Sy)
STAGGERED=np.kron(Sz,Sz)/2
UNIFORM=np.kron(I,Sz)/2

def hamiltonian(kx,ky):
    """Return (4x4 H, 4x4 dH/dkx); two exact PT doublets."""
    d=np.array([.92+.21*np.cos(kx)+.16*np.cos(ky),
                .24*np.sin(kx)+.17*np.sin(ky),
                .19*np.sin(kx)+.14*np.cos(ky),
                .18*np.sin(ky)+.15*np.cos(kx),
                .25+.19*np.cos(kx+ky)])
    dx=np.array([-.21*np.sin(kx),.24*np.cos(kx),.19*np.cos(kx),
                 -.15*np.sin(kx),-.19*np.sin(kx+ky)])
    d0=.13*np.cos(kx)+.08*np.sin(ky);d0x=-.13*np.sin(kx)
    return (d0*I4+np.einsum('a,aij->ij',d,G),
            d0x*I4+np.einsum('a,aij->ij',dx,G))

def projectors(kx,ky):
    H,_=hamiltonian(kx,ky)
    e,U=np.linalg.eigh(H)
    if max(abs(e[0]-e[1]),abs(e[2]-e[3]))>1.e-10:
        raise ValueError('model is not exactly two doublets')
    Pm=U[:,:2]@U[:,:2].conj().T
    Pp=I4-Pm
    return ((float((e[0]+e[1])/2),Pm),(float((e[2]+e[3])/2),Pp))

def projector_resolvent(kx,ky,X,eta=.1,order=1,omega=0.):
    """U(2)xU(2)-invariant PT-projector Liouvillian inverse."""
    if eta<=0 or order not in (1,2):raise ValueError('eta/order')
    X=np.asarray(X,dtype=complex)
    if X.shape!=(4,4):raise ValueError('X shape')
    out=np.zeros((4,4),complex)
    for ea,Pa in projectors(kx,ky):
        for eb,Pb in projectors(kx,ky):
            out+=Pa@X@Pb/(order*eta-1j*omega+1j*(ea-eb))
    return out

def transform_PT(X):
    return PT@np.asarray(X).conj()@PT.conj().T
