"""Conventional physical-spin-current, Cui et al. 2026 Appendix C (127)-(134).

H, v_a=dH/dk_a, K_ab=d²H/(dk_b dk_a) and S_c=physical Pauli spin
are in the SAME k-independent orthonormal basis. S_c is k-independent.
Natural units e=hbar=tau=1. O_a,c={v_a,S_c}/4 is conventional spin
current, not a generally conserved SOC spin current. RTA is model-specific.
No Wannier operator or material observable is inferred.
"""
from __future__ import annotations
import numpy as np
from .core import PAULI

def _validated(h, velocities, physical_spin, hessian=None, tol=1e-9):
    H=np.asarray(h,complex);v=np.asarray(velocities,complex);s=np.asarray(physical_spin,complex)
    if H.ndim != 2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError('H must be N x N, N>=2')
    n=H.shape[0]
    if v.shape!=(2,n,n) or s.shape!=(n,n):
        raise ValueError('Require v: (2,N,N), physical Pauli spin: (N,N)')
    items=[H,*v,s];K=None
    if hessian is not None:
        K=np.asarray(hessian,complex)
        if K.shape!=(2,2,n,n): raise ValueError('Require Hessian (2,2,N,N)')
        items.extend(K.reshape(-1,n,n))
        if not np.allclose(K,np.swapaxes(K,0,1),atol=tol,rtol=0):
            raise ValueError('Mixed momentum Hessian must be symmetric')
    for x in items:
        if not np.all(np.isfinite(x)) or not np.allclose(x,x.conj().T,atol=tol,rtol=0):
            raise ValueError('All operators must be finite and Hermitian')
    E,U=np.linalg.eigh(H)
    if np.min(np.diff(E))<=tol:
        raise ValueError('Degenerate bands require non-Abelian treatment')
    return H,v,s,K,E,U

def fermi(energies,mu=0.15,temperature=0.18):
    """f and df/dE at fixed mu. Zero-temperature derivative not provided."""
    if temperature<=0:raise ValueError('RTA derivative requires T>0')
    x=np.clip((np.asarray(energies)-mu)/temperature,-500,500)
    f=1/(1+np.exp(x))
    return f,-f*(1-f)/temperature

def spin_current_operator(velocities,spin):
    """Conventional spin-current operator {dH/dk_a,sigma_c}/4, a=x,y."""
    v=np.asarray(velocities,complex);s=np.asarray(spin,complex)
    return np.asarray([(va@s+s@va)/4 for va in v])

def intrinsic_spin_geometry(h,velocities,physical_spin,*,gap_tolerance=1e-9):
    """Per-band intrinsic spin Kubo plus Eq.(128) gA/OmegaN sector oracle.

    C_nm^ab = O_a,nm v_b,mn / (E_n-E_m)^2.
    Intrinsic = -2 Im sum_m C; sectors give -(2*gA - OmegaN).
    """
    _,v,s,_,E,U=_validated(h,velocities,physical_spin,tol=gap_tolerance)
    O=spin_current_operator(v,s)
    vb=np.einsum('in,aij,jm->anm',U.conj(),v,U,optimize=True)
    ob=np.einsum('in,aij,jm->anm',U.conj(),O,U,optimize=True)
    N=len(E);T=np.zeros((N,2,2),complex)
    for n in range(N):
        for m in range(N):
            if n==m:continue
            T[n]+=np.outer(ob[:,n,m],vb[:,m,n])/(E[n]-E[m])**2
    gA=((T+np.swapaxes(T,1,2))/2).imag
    omegaN=-(T-np.swapaxes(T,1,2)).imag
    return dict(energies=E,eigvecs=U,velocity_bands=vb,spin_current_bands=ob,
                complex_spin_qgt=T,gA=gA,omegaN=omegaN,
                intrinsic=-2*T.imag,sector_intrinsic=-(2*gA-omegaN))

def spin_current_derivative(h,velocities,hessian,physical_spin,*,gap_tolerance=1e-9):
    """Exact -partial_b <O_a> for fixed physical spin, including Hessian.

    partial_b<O_a>_n = 2 Re sum_m O_a,nm v_b,mn/(E_n-E_m)
                      + <n|{d_b d_a H,S_c}/4|n>.
    Geometric-only Eq.(133) needs vanishing explicit-vertex expectation.
    """
    _,v,s,K,E,U=_validated(h,velocities,physical_spin,hessian,gap_tolerance)
    if K is None: raise ValueError('Full derivative requires explicit Hessian')
    q=intrinsic_spin_geometry(h,v,s,gap_tolerance=gap_tolerance)
    vb=q['velocity_bands'];ob=q['spin_current_bands'];N=len(E)
    geometric=np.zeros((N,2,2),float)
    for n in range(N):
        for m in range(N):
            if n==m:continue
            geometric[n]-=2*np.outer(ob[:,n,m],vb[:,m,n]).real/(E[n]-E[m])
    vertex=spin_current_operator(K.reshape(4,N,N),s).reshape((2,2,N,N))
    vertex_exp=np.array([[[np.vdot(U[:,n],vertex[a,b]@U[:,n]).real
                            for b in range(2)] for a in range(2)] for n in range(N)])
    return {**q,'geometry_only':geometric,'hessian_vertex':vertex_exp,
            'full_RTA':geometric-vertex_exp,
            'band_spin_current':np.stack([np.diag(ob[a]).real for a in range(2)],axis=1),
            'band_velocity':np.stack([np.diag(vb[b]).real for b in range(2)],axis=1)}

def spin_orbital_benchmark(kx,ky,*,periodic=True):
    """Synthetic four-band orbital-tau x physical-spin-sigma benchmark."""
    I=np.eye(2,dtype=complex);x,y,z=PAULI;K=lambda a,b:np.kron(a,b)
    opx=K(x,I)+.19*K(y,x);opy=K(y,I)+.16*K(x,y);mass=K(z,I)
    H0=.7*mass+.27*K(z,x)+.36*K(I,z)
    if periodic:
        H=(H0+np.sin(kx)*opx+np.sin(ky)*opy
           +.35*(np.cos(kx)+np.cos(ky)-2)*mass)
        V=np.stack([np.cos(kx)*opx-.35*np.sin(kx)*mass,
                    np.cos(ky)*opy-.35*np.sin(ky)*mass])
        D2=np.zeros((2,2,4,4),complex)
        D2[0,0]=-np.sin(kx)*opx-.35*np.cos(kx)*mass
        D2[1,1]=-np.sin(ky)*opy-.35*np.cos(ky)*mass
    else:
        H=H0+kx*opx+ky*opy
        V=np.stack([opx,opy]);D2=np.zeros((2,2,4,4),complex)
    return H,V,D2,K(I,z)

def integrate_spin_rta_bz(*,mesh=40,mu=.15,temperature=.18,spin=None):
    """BZ midpoint, tau=1: direct f'<O>v vs exact -f d<O> and reduced -f geometric.

    Outputs dimensionless 2D BZ integral normalized by (2pi)^2.
    """
    if mesh<4:raise ValueError('mesh must be >=4')
    accum=np.zeros((3,2,2),float)
    for kx in -np.pi+(np.arange(mesh)+.5)*2*np.pi/mesh:
        for ky in -np.pi+(np.arange(mesh)+.5)*2*np.pi/mesh:
            h,v,K,s=spin_orbital_benchmark(kx,ky,periodic=True)
            if spin is not None:s=spin
            d=spin_current_derivative(h,v,K,s)
            f,fp=fermi(d['energies'],mu,temperature)
            accum[0]+=np.einsum('n,na,nb->ab',fp,d['band_spin_current'],d['band_velocity'])
            accum[1]+=np.einsum('n,nab->ab',f,d['full_RTA'])
            accum[2]+=np.einsum('n,nab->ab',f,d['geometry_only'])
    direct,full,geom=accum/mesh**2
    return dict(mesh=mesh,direct_fermi_surface=direct,full_fermi_sea=full,
                geometry_only=geom,missing_hessian=full-geom,
                normalization='2D BZ d2k/(2pi)^2; tau=1; conventional spin current')
