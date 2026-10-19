"""Reference Kane–Mele and small test Hamiltonians.

The honeycomb convention reproduces the independently audited sign choice of
Round 16; this is a convention, NOT an inferred author raw-data setting.
"""
from __future__ import annotations
import numpy as np
from .core import BlochModel
I2=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1,-1]).astype(complex)


def kane_mele(*,bond_nm=.142, t=.85, intrinsic_soc=.010,rashba_soc=.020,
              staggered_potential=.010,g=10.):
    """Graphene-style two-sublattice, spinful 4x4 Bloch model.

    sigma_z sign and Rashba orientation calibrated against the previous
    independent Kane–Mele test, NOT claimed as authors' unique convention.
    """
    nn=np.array([[0,1],[np.sqrt(3)/2,-.5],[-np.sqrt(3)/2,-.5]])*bond_nm
    a=np.stack((nn[0]-nn[1],nn[0]-nn[2]))
    G=2*np.pi*np.linalg.inv(a).T
    nnn=np.array([nn[1]-nn[2],nn[2]-nn[0],nn[0]-nn[1]])
    # With our choice of directed bonds, paper +SOC = -SOC in legacy H.
    lso=-intrinsic_soc; lr=-rashba_soc; lv=staggered_potential
    spin={'x':np.kron(I2,SX),'y':np.kron(I2,SY),'z':np.kron(I2,SZ)}
    def blocks(k, axis=None):
        k=np.atleast_2d(np.asarray(k,float));nk=len(k)
        ek=np.exp(1j*k@nn.T);sok=k@nnn.T
        so=2*lso*np.sin(sok).sum(axis=1)
        # Offdiagonal from A->B bond harmonics.
        hop=-t*ek.sum(1)[:,None,None]*I2
        if axis is None:
            vh=None
        else:
            vh=-1j*t*np.sum(ek*nn[None,:,axis],axis=1)[:,None,None]*I2
        for j in range(3):
            orient=SX*nn[j,1]/bond_nm-SY*nn[j,0]/bond_nm
            hop+=1j*lr*ek[:,j,None,None]*orient
            if axis is not None:vh+=-lr*nn[j,axis]*ek[:,j,None,None]*orient
        h=np.zeros((nk,4,4),complex)
        if axis is None:
            zz=so
        else:
            zz=2*lso*np.sum(np.cos(sok)*nnn[None,:,axis],axis=1)
        h[:,:2,:2]=zz[:,None,None]*SZ
        h[:,2:,2:]=-zz[:,None,None]*SZ
        if axis is None:
            h[:,:2,:2]+=lv*I2
            h[:,2:,2:]-=lv*I2
        h[:,:2,2:]=hop if axis is None else vh
        h[:,2:,:2]=h[:,:2,2:].conj().swapaxes(1,2)
        return h
    ops={name:lambda k,s=s:np.broadcast_to(s,(len(np.atleast_2d(k)),4,4)).copy() for name,s in spin.items()}
    return BlochModel(lambda k:blocks(k),lambda k,beta:blocks(k,beta),ops,G,g_factor=g,name='KM-4-band')


def dirac_1d(*,m=.3,delta=.5):
    """1D two-level analytic benchmark; k is dimensionless test coordinate."""
    G=np.array([[2*np.pi]])
    def h(k):
        x=np.asarray(k)[:,0];return ((m+np.cos(x))[:,None,None]*SX+
                                     np.sin(x)[:,None,None]*SY+delta*SZ)
    def velocity(k,beta):
        x=np.asarray(k)[:,0];return -np.sin(x)[:,None,None]*SX+np.cos(x)[:,None,None]*SY
    ops={a:lambda k,s=s:np.broadcast_to(s,(len(k),2,2)).copy() for a,s in {'x':SX,'y':SY,'z':SZ}.items()}
    return BlochModel(h,velocity,ops,G,name='dirac-1D-analytic')
