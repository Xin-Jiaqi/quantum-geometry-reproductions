"""Gauge-covariant interband response for exact/near-degenerate band clusters.

This is a cluster-flattened clean interband protocol, not a conserving
metallic Kubo formula. Degenerate states are treated as a projector/subspace;
intra-cluster position matrix elements are set to zero because the chosen
interband-only response excludes them. A cluster within gap_tolerance is
flattened to its mean energy before both Kubo and Liouville evaluations.
"""
from __future__ import annotations
import numpy as np
from .causal_response import _fermi


def _cluster_data(H,dH,spin,*,gap_tolerance=1e-9):
    H=np.asarray(H,complex);dH=np.asarray(dH,complex);spin=np.asarray(spin,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError("H must be square N>=2")
    N=H.shape[0]
    if dH.shape!=(2,N,N) or spin.shape!=(3,N,N):
        raise ValueError("Require dH=(2,N,N), physical spin=(3,N,N)")
    if not np.isfinite(gap_tolerance) or gap_tolerance<0:
        raise ValueError("gap_tolerance must be finite and nonnegative")
    for A in (H,*dH,*spin):
        if not np.all(np.isfinite(A)):raise ValueError("Nonfinite input")
        if not np.allclose(A,A.conj().T,rtol=0,atol=1e-10):
            raise ValueError("All operators must be Hermitian")
    e,U=np.linalg.eigh(H)
    scale=max(1.,float(np.max(abs(e))))
    tol=gap_tolerance*scale
    labels=np.zeros(N,dtype=int)
    for n in range(1,N):
        labels[n]=labels[n-1]+int(e[n]-e[n-1]>tol)
    centers=np.array([np.mean(e[labels==g]) for g in range(labels[-1]+1)])
    ec=centers[labels]
    v=np.einsum("ni,anm,mj->aij",U.conj(),dH,U,optimize=True)
    s=np.einsum("ni,bnm,mj->bij",U.conj(),spin,U,optimize=True)
    r=np.zeros_like(v)
    for n in range(N):
        for m in range(N):
            if labels[n]!=labels[m]:
                r[:,n,m]=v[:,n,m]/(1j*(ec[n]-ec[m]))
    if not np.allclose(r,np.swapaxes(r.conj(),-1,-2),rtol=0,atol=2e-10):
        raise AssertionError("Inter-cluster position operator is not Hermitian")
    R=np.einsum("in,anm,jm->aij",U,r,U.conj(),optimize=True)
    Heff=(U*ec)@U.conj().T
    return dict(e=ec,raw_e=e,U=U,labels=labels,centers=centers,
                v=v,s=s,r=r,R=R,Heff=Heff,scale=scale,tolerance=tol)


def cluster_structure(H,*,gap_tolerance=1e-9):
    """Return energy centers and basis-invariant spectral projectors per cluster."""
    H=np.asarray(H,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError("H must be square N>=2")
    e,U=np.linalg.eigh(H);scale=max(1.,float(np.max(abs(e))))
    if gap_tolerance<0 or not np.isfinite(gap_tolerance):
        raise ValueError("gap_tolerance must be finite and nonnegative")
    labels=np.zeros(len(e),dtype=int);tol=gap_tolerance*scale
    for n in range(1,len(e)):labels[n]=labels[n-1]+int(e[n]-e[n-1]>tol)
    centers=np.array([np.mean(e[labels==g]) for g in range(labels[-1]+1)])
    projectors=[]
    for g in range(len(centers)):
        W=U[:,labels==g];projectors.append(W@W.conj().T)
    return dict(centers=centers,labels=labels,projectors=projectors,
                dimensions=[int(np.sum(labels==g)) for g in range(len(centers))])


def cluster_zeeman_qgt(H,dH,spin,*,cluster_index=0,gap_tolerance=1e-9):
    r"""Non-Abelian Zeeman QGT block T_C^{ab}=P_C r^a (1-P_C) s^b P_C.

    Returned matrices are expressed in the numerically chosen basis of cluster
    C. Under a change of basis W_C they transform as W_C^\dag T W_C; traces
    and eigenvalues are gauge invariant. This is a geometric block, not SI.
    """
    D=_cluster_data(H,dH,spin,gap_tolerance=gap_tolerance)
    groups=len(D["centers"])
    if not isinstance(cluster_index,(int,np.integer)) or not 0<=cluster_index<groups:
        raise ValueError("cluster_index out of range")
    idx=np.flatnonzero(D["labels"]==cluster_index)
    ext=np.flatnonzero(D["labels"]!=cluster_index)
    T=np.zeros((2,3,len(idx),len(idx)),complex)
    if len(ext):
        for a in range(2):
            for b in range(3):
                T[a,b]=D["r"][a][np.ix_(idx,ext)]@D["s"][b][np.ix_(ext,idx)]
    return dict(energy=float(D["centers"][cluster_index]),dimension=len(idx),
                qgt=T,trace=np.trace(T,axis1=-2,axis2=-1),
                projector=cluster_structure(H,gap_tolerance=gap_tolerance)["projectors"][cluster_index])


def cluster_interband_kubo(H,dH,spin,*,mu=0.,temperature=.025,z=.04,kappa=1.,
                           gap_tolerance=1e-9):
    """Gauge-invariant directed cluster-pair Kubo response (2x3, 3x2)."""
    D=_cluster_data(H,dH,spin,gap_tolerance=gap_tolerance)
    e,v,s,r,labels=(D[k] for k in ("e","v","s","r","labels"))
    f=_fermi(e,mu,temperature);z=complex(z)
    JB=np.zeros((2,3),complex);ME=np.zeros((3,2),complex)
    for n in range(len(e)):
        for m in range(len(e)):
            if labels[n]==labels[m]:continue
            q=f[n]-f[m];den=z-(e[m]-e[n])
            if abs(den)<1e-10:
                if abs(q)>1e-10:raise ValueError("On-shell pole")
                continue
            JB-=kappa*q*np.outer(v[:,n,m],s[:,m,n])/den
            ME+=kappa*q*np.outer(s[:,n,m],r[:,m,n])/den
    return JB,ME


def cluster_interband_liouville(H,dH,spin,*,mu=0.,temperature=.025,z=.04,kappa=1.,
                                gap_tolerance=1e-9):
    """Independent orbital-basis N² Liouville solve for the cluster-flat model."""
    D=_cluster_data(H,dH,spin,gap_tolerance=gap_tolerance)
    Hc=D["Heff"];U=D["U"];e=D["e"];R=D["R"];z=complex(z);N=len(e)
    if abs(z)<1e-12:raise ValueError("Full Liouvillian singular at z=0")
    f=_fermi(e,mu,temperature);rho0=(U*f)@U.conj().T
    L=np.kron(np.eye(N),Hc)-np.kron(Hc.T,np.eye(N))
    M=z*np.eye(N*N,dtype=complex)-L
    dH=np.asarray(dH,complex);spin=np.asarray(spin,complex)
    JB=np.zeros((2,3),complex);ME=np.zeros((3,2),complex)
    for b,S in enumerate(spin):
        Q=-kappa*(S@rho0-rho0@S)
        d=np.linalg.solve(M,Q.reshape(-1,order="F")).reshape((N,N),order="F")
        for a,V in enumerate(dH):JB[a,b]=np.trace(V@d)
    for a,P in enumerate(R):
        Q=P@rho0-rho0@P
        d=np.linalg.solve(M,Q.reshape(-1,order="F")).reshape((N,N),order="F")
        for b,S in enumerate(spin):ME[b,a]=kappa*np.trace(S@d)
    return JB,ME
