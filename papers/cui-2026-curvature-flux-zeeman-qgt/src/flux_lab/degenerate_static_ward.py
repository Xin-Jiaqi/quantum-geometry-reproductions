"""Round24: exact finite-T static Zeeman-current Ward identity with degeneracies.

H_B(k)=H(k)-kappa*B_b*S_b and J_a=dH/dk_a with k-independent S_b.
At finite temperature the Fréchet derivative of Fermi occupation is

    Df_H[A] = sum_nm f[eps_n,eps_m] |n><n|A|m><m|,

with f[e,e]=f'(e), INCLUDING n!=m within an exactly degenerate cluster.
Therefore chi_JB^ab=-kappa Tr[V_a Df_H[S_b]]
                  =-kappa d_ka Tr[f(H) S_b].

No arbitrary choice of diagonal eigenvectors inside Kramers doublets enters
the complete response. Cluster decomposition is a diagnostic for exact
degeneracies only; it does NOT flatten near-degenerate energies or approximate
the response. This is STATIC EQUILIBRIUM, not a finite-frequency metallic KME.
"""
from __future__ import annotations
import numpy as np
from .core import PAULI
from .causal_response import _fermi
from .hexagonal_bz import BASIS,INV_BASIS


def _input(H,V,S,mu,temperature,kappa):
    H=np.asarray(H,complex);V=np.asarray(V,complex);S=np.asarray(S,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError("H must be square N>=2")
    n=len(H)
    if V.shape!=(2,n,n) or S.shape!=(3,n,n):
        raise ValueError("V=(2,N,N) and physical spin S=(3,N,N) required")
    if not (np.isfinite(mu) and np.isfinite(temperature) and temperature>0
            and np.isfinite(kappa)):
        raise ValueError("finite mu, kappa and strictly positive T required")
    for A in (H,*V,*S):
        if not np.all(np.isfinite(A)) or not np.allclose(A,A.conj().T,rtol=0,atol=1e-10):
            raise ValueError("finite Hermitian H, velocity, spin required")
    return H,V,S


def fermi_divided_difference(energies,*,mu=0.,temperature=.2):
    """Stable symmetric f[e_i,e_j], with exact repeated-energy limit f'(e).

    The smooth small-gap hyperbolic formula handles near degeneracies; at
    larger gaps ordinary Fermi divided differences avoid cosh overflow.
    """
    e=np.asarray(energies,float)
    if e.ndim!=1 or not np.all(np.isfinite(e)) or e.size<2:
        raise ValueError("finite real energies (N>=2) required")
    if not np.isfinite(temperature) or temperature<=0 or not np.isfinite(mu):
        raise ValueError("finite mu and strictly positive T required")
    f=_fermi(e,mu,temperature)
    de=e[:,None]-e[None,:]
    close=np.abs(de)<.1*temperature
    out=np.empty_like(de)
    np.divide(f[:,None]-f[None,:],de,out=out,where=~close)
    x=np.clip((e-mu)/(2*temperature),-300,300)
    d=de/(2*temperature)
    # sinc_h(d)=sinh(d)/d is analytic, including d=0.
    ratio=np.ones_like(de)
    sel=close & (np.abs(d)>1e-9)
    ratio[sel]=np.sinh(d[sel])/d[sel]
    denom=4*temperature*np.cosh(x)[:,None]*np.cosh(x)[None,:]
    out[close]=-ratio[close]/denom[close]
    return out


def degenerate_static_ward(H,V,S,*,mu=0.,temperature=.2,kappa=1.,
                           degeneracy_tolerance=1e-12):
    """Full static chi and a basis-stable exact-degeneracy sector partition.

    For fixed eigenvalues the f-divided-difference matrix is regular even
    when every band is degenerate. The total is gauge covariant. Diagnostic
    intra/inter flags use a declared absolute+relative energy tolerance and
    are canonical ONLY when the subspaces are exactly degenerate.
    """
    H,V,S=_input(H,V,S,mu,temperature,kappa)
    if not np.isfinite(degeneracy_tolerance) or degeneracy_tolerance<=0:
        raise ValueError("positive finite degeneracy tolerance required")
    e,U=np.linalg.eigh(H)
    v=np.einsum("mi,amn,nj->aij",U.conj(),V,U,optimize=True)
    s=np.einsum("mi,bmn,nj->bij",U.conj(),S,U,optimize=True)
    dd=fermi_divided_difference(e,mu=mu,temperature=temperature)
    scale=max(1.,float(np.max(np.abs(e))))
    mask=np.abs(e[:,None]-e[None,:]) <= degeneracy_tolerance*scale
    # v[a,m,n] s[b,n,m] f[eps_m,eps_n]
    contract=lambda kernel: -kappa*np.einsum(
        "amn,bnm,mn->ab",v,s,kernel,optimize=True)
    intra=contract(dd*mask)
    inter=contract(dd*(~mask))
    full=contract(dd)
    if not np.allclose(intra+inter,full,rtol=2e-12,atol=2e-12):
        raise AssertionError("divided-difference partition inconsistent")
    f=_fermi(e,mu,temperature)
    rho=(U*f)@U.conj().T
    spin_density=np.array([np.trace(rho@x) for x in S])
    dims=[];start=0
    for i in range(1,len(e)+1):
        if i==len(e) or not mask[i-1,i]:
            dims.append(i-start);start=i
    return dict(full_static=full,intracluster=intra,intercluster=inter,
                energies=e,cluster_dimensions=dims,
                spin_density=spin_density,fermi_divided_difference=dd,
                status="complete finite-temperature equilibrium static response")


def projector_static_ward(H,V,S,*,mu=0.,temperature=.2,kappa=1.,
                          degeneracy_tolerance=1e-12):
    """Independent spectral-projector contraction (no band-index matrix sum).

    Only EXACT degeneracies admit an unambiguous intra-cluster projector
    sector. For numerical eigenspectrum differences within declared tol,
    centers are used in this ORACLE, so it should not be used for broad
    near-degenerate flattening as a physical approximation.
    """
    H,V,S=_input(H,V,S,mu,temperature,kappa)
    e,U=np.linalg.eigh(H);scale=max(1.,np.max(abs(e)))
    groups=[];start=0
    for i in range(1,len(e)+1):
        if i==len(e) or e[i]-e[i-1]>degeneracy_tolerance*scale:
            W=U[:,start:i]
            groups.append((float(np.mean(e[start:i])),W@W.conj().T))
            start=i
    centers=np.array([a for a,_ in groups])
    divided=fermi_divided_difference(centers,mu=mu,temperature=temperature)
    out=np.zeros((2,3),complex)
    intra=np.zeros((2,3),complex)
    for i,(ei,Pi) in enumerate(groups):
        for j,(ej,Pj) in enumerate(groups):
            for a in range(2):
                for b in range(3):
                    item=-kappa*divided[i,j]*np.trace(Pi@V[a]@Pj@S[b])
                    out[a,b]+=item
                    if i==j:intra[a,b]+=item
    return dict(full_static=out,intracluster=intra,
                intercluster=out-intra,cluster_dimensions=[int(round(np.trace(P).real)) for _,P in groups])


def kramers_periodic_model(k):
    """Triangular-periodic four-band TR-symmetric spin-orbital test.

    Kramers degeneracy is exact at k=0, where every velocity matrix can
    couple states inside an exactly degenerate doublet. S is genuine
    I_tau x sigma and independent of k. Derivatives returned in Cartesian k.
    """
    k=np.asarray(k,float)
    if k.shape!=(2,) or not np.all(np.isfinite(k)):
        raise ValueError("k must be finite Cartesian vector length 2")
    q=INV_BASIS@k
    t=PAULI
    Id=np.eye(2,dtype=complex)
    tau=lambda A,B:np.kron(A,B)
    Z=tau(t[2],Id)
    W1=.55*tau(t[0],t[0])+.21*tau(t[1],Id)
    W2=.47*tau(t[0],t[1])+.13*tau(t[0],t[2])
    H=(1.2*Z+.31*tau(t[1],t[0])
       +np.sin(q[0])*W1+np.sin(q[1])*W2
       +.12*(np.cos(q[0])+np.cos(q[1])-2)*Z)
    Vq=np.stack([np.cos(q[0])*W1-.12*np.sin(q[0])*Z,
                 np.cos(q[1])*W2-.12*np.sin(q[1])*Z])
    V=np.einsum("ab,bij->aij",INV_BASIS.T,Vq,optimize=True)
    S=np.stack([tau(Id,s) for s in t])
    return H,V,S


def finite_zeeman_static(H,V,S,*,mu=0.,temperature=.2,kappa=1.,step=1e-5):
    """Finite uniform B derivative; no divided-difference tensor used."""
    from .equilibrium_ward import equilibrium_current
    H,V,S=_input(H,V,S,mu,temperature,kappa)
    if step<=0 or not np.isfinite(step):raise ValueError("positive step required")
    out=np.zeros((2,3),complex)
    for b in range(3):
        plus=equilibrium_current(H-kappa*step*S[b],V,mu=mu,temperature=temperature)
        minus=equilibrium_current(H+kappa*step*S[b],V,mu=mu,temperature=temperature)
        out[:,b]=(plus-minus)/(2*step)
    return out


def finite_k_spin_gradient(k,*,mu=0.,temperature=.2,kappa=1.,step=1e-5):
    """Finite k derivative of physical spin density for periodic toy model."""
    from .equilibrium_ward import spin_density
    k=np.asarray(k,float)
    if k.shape!=(2,) or step<=0 or not np.isfinite(step):
        raise ValueError("finite k length 2 and positive step required")
    out=np.zeros((2,3),complex)
    for a in range(2):
        dk=np.zeros(2);dk[a]=step
        hp,_,sp=kramers_periodic_model(k+dk)
        hm,_,sm=kramers_periodic_model(k-dk)
        out[a]=-kappa*(spin_density(hp,sp,mu=mu,temperature=temperature)
                       -spin_density(hm,sm,mu=mu,temperature=temperature))/(2*step)
    return out
