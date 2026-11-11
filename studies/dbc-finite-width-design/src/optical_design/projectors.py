"""Independent projector-derivative reference, not a production Wannier engine.

Cproj^{mn}_{a;bc}=Tr[P_n d_b P_m (d_ac P_n + d_a P_m d_c P_n)]
from Mitscherling et al., PRB 112, 085104 (2025), DOI qscv-qxqt.
For rank-one bands Cproj = - C_Zou with the n,m indexing used here.
Uses fourth-order finite differences of projectors, without phase alignment.
Exactly degenerate, isolated groups are allowed; split groups are never flattened.
"""
from __future__ import annotations
import numpy as np


def projector_kernel(hamiltonian, k, nocc=1, axes=(0,0,1), *, groups=None,
                     step=1e-3, gap_tolerance=1e-9):
    k=np.asarray(k,float)
    if k.ndim!=1 or not np.isfinite(k).all() or not np.isfinite(step) or step<=0:
        raise ValueError('Finite vector k and positive step required')
    if not np.isfinite(gap_tolerance) or gap_tolerance<=0:
        raise ValueError('Positive gap tolerance required')
    ndim=len(k);a,b,c=axes
    if any(type(x) is not int or x<0 or x>=ndim for x in axes):
        raise ValueError('Axes must index the momentum coordinates')
    H=np.asarray(hamiltonian(k),complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1]:raise ValueError('Square Hamiltonian required')
    N=len(H)
    if type(nocc) is not int or not 0<nocc<N:raise ValueError('Invalid occupation')
    groups=[(i,) for i in range(N)] if groups is None else [tuple(g) for g in groups]
    flat=[i for g in groups for i in g]
    if sorted(flat)!=list(range(N)) or any(not g for g in groups):
        raise ValueError('Groups must partition all bands exactly once')
    if any(min(g)<nocc<=max(g) for g in groups):raise ValueError('A group straddles the occupation gap')
    cache={}
    def sample(offset):
        key=tuple(offset)
        if key not in cache:
            h=np.asarray(hamiltonian(k+step*np.asarray(key)),complex)
            if h.shape!=(N,N) or not np.isfinite(h).all() or not np.allclose(h,h.conj().T,atol=1e-12,rtol=1e-12):
                raise ValueError('Hamiltonian must stay finite and Hermitian')
            e,u=np.linalg.eigh(h);ec=[];P=[]
            for g in groups:
                eg=e[list(g)]
                if np.ptp(eg)>gap_tolerance:
                    raise ValueError('Split group: flattening would change the physical problem')
                ec.append(float(np.mean(eg)));v=u[:,list(g)];P.append(v@v.conj().T)
            ec=np.asarray(ec)
            for i in range(len(ec)):
                if np.any(abs(ec[i+1:]-ec[i])<gap_tolerance):
                    raise ValueError('Unresolved inter-group degeneracy')
            cache[key]=(ec,np.asarray(P))
        return cache[key]
    zero=(0,)*ndim;energies,P=sample(zero)
    offsets=(-2,-1,1,2);weights=np.array([1,-8,8,-1])/(12*step)
    dP=np.zeros((ndim,)+P.shape,complex)
    ddP=np.zeros((ndim,ndim)+P.shape,complex)
    for i in range(ndim):
        for t,w in zip(offsets,weights):
            delta=np.zeros(ndim,int);delta[i]=t;dP[i]+=w*sample(delta)[1]
        for t,w in zip((-2,-1,0,1,2),(-1,16,-30,16,-1)):
            delta=np.zeros(ndim,int);delta[i]=t;ddP[i,i]+=w*sample(delta)[1]/(12*step**2)
        for j in range(i):
            for t,w in zip(offsets,weights):
                for u,v in zip(offsets,weights):
                    delta=np.zeros(ndim,int);delta[i]=t;delta[j]=u
                    ddP[i,j]+=w*v*sample(delta)[1]
            ddP[j,i]=ddP[i,j]
    gaps=[];kernels=[];pairs=[]
    for n,gn in enumerate(groups):
        if min(gn)>=nocc:continue
        for m,gm in enumerate(groups):
            if max(gm)<nocc:continue
            def C(b,c):return np.trace(P[n]@dP[b,m]@(ddP[a,c,n]+dP[a,m]@dP[c,n]))
            kernels.append(float(np.imag(C(b,c)+C(c,b))))
            gaps.append(float(energies[m]-energies[n]));pairs.append((n,m))
    return {'gaps':np.asarray(gaps),'kernels':np.asarray(kernels),'pairs':pairs,
            'groups':groups,'step':step,'sample_count':len(cache),
            'scope':'Fixed orthonormal basis, isolated scalar-energy groups, T=0 insulator.'}
