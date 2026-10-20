"""Gauge-local valley k.p expansion and occupied-to-empty STATIC response.

Taylor coefficients are constructed from the identical 8-band microscopic model;
no fitted response or presumed Dirac mass. The orbital embedding is inherited
from its Bloch gauge, and V_kp = dH_kp/dq is enforced algebraically.
A 2D analytic Taylor patch is NOT a globally periodic material Hamiltonian.
"""
from dataclasses import dataclass
import numpy as np
from .sliding_tight_binding import SlidingBilayer
from nmee_lab.clean_core.core import MU_B

@dataclass(frozen=True)
class ValleyKP:
    center:np.ndarray
    h0:np.ndarray
    v0:np.ndarray  # (2,8,8)
    h2:np.ndarray  # (2,2,8,8)

    @classmethod
    def from_microscopic(cls, model: SlidingBilayer, center, fd=1e-3):
        center=np.asarray(center,float)
        if center.shape!=(2,) or not fd>0:raise ValueError('bad k.p origin/derivative step')
        h=model.matrices(center[None,:])[0]
        v=np.stack([model.matrices(center[None,:],axis=j)[0] for j in range(2)])
        h2=np.zeros((2,2,len(h),len(h)),complex)
        for i in range(2):
            displacement=np.eye(2)[i]*fd
            for j in range(2):
                plus=model.matrices((center+displacement)[None,:],axis=j)[0]
                minus=model.matrices((center-displacement)[None,:],axis=j)[0]
                h2[i,j]=(plus-minus)/(2*fd)
        h2=(h2+h2.swapaxes(0,1))/2
        return cls(center,h,v,h2)

    def matrices(self,k,axis=None,order=2):
        k=np.atleast_2d(np.asarray(k,float))
        q=k-self.center
        if order not in (1,2):raise ValueError('Taylor order must be 1 or 2')
        if axis is None:
            H=self.h0[None,:,:]+np.einsum('ni,ijk->njk',q,self.v0,optimize=True)
            if order==2:H=H+.5*np.einsum('ni,nj,ijab->nab',q,q,self.h2,optimize=True)
            return H
        if axis not in (0,1):raise ValueError('axis must be 0/1')
        V=np.broadcast_to(self.v0[axis],(len(k),*self.v0[axis].shape)).copy()
        if order==2:V+=np.einsum('ni,iab->nab',q,self.h2[:,axis],optimize=True)
        return V


def static_kernel(H,V,spins, sz, mu, *, g_factor=10, step_T=.25):
    """Integrand (N,2,2) for d_B of occupied spin-electric mixed curvature."""
    H=np.asarray(H);V=np.asarray(V)
    if H.ndim!=3 or H.shape[1:]!=(8,8) or V.shape!=(2,len(H),8,8):
        raise ValueError('incorrect H or dH/dk shapes')
    e0=np.linalg.eigvalsh(H)
    if np.any(e0[:,3]>=mu) or np.any(e0[:,4]<=mu):
        raise ValueError('insulating occupation gap closed')
    out=[]
    for sg in (1,-1):
        es,U=np.linalg.eigh(H-sg*g_factor*MU_B*step_T*sz)
        if np.any(es[:,3]>=mu) or np.any(es[:,4]<=mu):raise ValueError('field closes gap')
        Ud=U.conj().swapaxes(1,2)
        s=[(Ud@S@U)[:,:4,4:] for S in spins]
        vs=[np.swapaxes((Ud@v@U)[:,4:,:4],1,2) for v in V]
        den=(es[:,:4,None]-es[:,None,4:])**2
        out.append(np.stack([np.stack([-2*np.imag(S*vel/den).sum(axis=(1,2)) for vel in vs],axis=1)
                             for S in s],axis=1))
    return (out[0]-out[1])/(2*step_T)


def compare_patch(spec,k,weights,centers,point_valley,*,order=2,step_T=.25,batch=1024,
                  same_scale=1.,cross_scale=1.):
    """Returns full and k.p valley-integrated 2x2 NMEE response."""
    obj=SlidingBilayer(spec)
    obj.hops=[(a,b,d,amp*(same_scale if a==b else cross_scale)) for a,b,d,amp in obj.hops]
    kp=[ValleyKP.from_microscopic(obj,c) for c in centers]
    spins=[obj.spins['total_x'],obj.spins['total_y']]
    sz=obj.spins['total_z']
    ans=np.zeros((2,2,2),float) # full vs k.p
    mine=np.full(2,np.inf)
    for v in range(2):
        kv=k[v*point_valley:(v+1)*point_valley]
        wv=weights[v*point_valley:(v+1)*point_valley]
        for i in range(0,len(kv),batch):
            sub=kv[i:i+batch]; w=wv[i:i+batch]
            for kind in range(2):
                H=obj.matrices(sub) if kind==0 else kp[v].matrices(sub,order=order)
                V=np.stack([obj.matrices(sub,axis=j) if kind==0 else kp[v].matrices(sub,axis=j,order=order) for j in range(2)])
                ker=static_kernel(H,V,spins,sz,spec.mu_eV,g_factor=obj.mono.g_factor,step_T=step_T)
                ans[kind]+=np.einsum('nij,n->ij',ker,w)
                e=np.linalg.eigvalsh(H)
                mine[kind]=min(mine[kind],float(np.min(e[:,4]-spec.mu_eV)),float(np.min(spec.mu_eV-e[:,3])))
    return {'full':ans[0], 'kp':ans[1], 'diff':ans[1]-ans[0], 'min_mu_gap':mine}
