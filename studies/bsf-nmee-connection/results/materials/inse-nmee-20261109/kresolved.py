"""Uniform-weight k-resolved form of the existing full-KS spin response."""
import numpy as np
from analyze import adj

def c_k(r,h,b=.002):
    pieces=[]
    for sign in (1,-1):
        en,u=np.linalg.eigh(h+sign*b*r.z);occ=u[:,:,:r.nocc];emp=u[:,:,r.nocc:]
        deriv=np.zeros_like(occ,dtype=complex)
        for a in range(4):
            v=r.links[:,a]@occ[r.neighbors[:,a]];o=adj(occ)@v
            deriv+=(-1)**a*r.mesh*r.cell[a//2,0]/(4*np.pi)*v@np.linalg.inv(o)
        dip=adj(emp)@(1j*deriv);delta=en[:,:r.nocc,None]-en[:,None,r.nocc:]
        total,lower=[(2*((adj(occ)@op@emp)*dip.swapaxes(1,2)).real/delta).sum(axis=(1,2))*r.factor for op in r.ops]
        pieces.append(np.stack([total,lower,total-lower],axis=1))
    return (pieces[0]-pieces[1])/(2*b)
