"""Round23: equilibrium static Zeeman-current Ward/no-current identity.

For fixed k-independent physical spin S_b and a smooth periodic Bloch H(k),
H_B(k) = H(k) - kappa B_b S_b, J_a = dH/dk_a, rho=f(H_B).

The static response including interband AND Fermi-occupation derivative is

 chi_{JB}^{ab}(k,0) =
 -kappa sum_{n!=m} f[e_n,e_m] v^a_nm S^b_mn
 -kappa sum_n f'(e_n) v^a_nn S^b_nn
 = -kappa partial_{k_a} Tr[f(H(k)) S_b].

Here f[e_n,e_m]=(f_n-f_m)/(e_n-e_m), and f'=df/de.
For a periodic BZ, integral of this derivative is zero (equilibrium current).
The individual interband term need NOT integrate to zero in a metal; the
Fermi-surface occupation-derivative term is essential.

These are synthetic clean equilibrium Bloch identities. This module does
NOT supply finite-frequency metallic KME, collisions, transverse contacts,
a nontrivial k-dependent Wannier spin embedding, or SI normalization.
"""
from __future__ import annotations
import numpy as np
from .causal_response import _fermi


def static_zeeman_current_ward(H,V,S,*,mu=0.,temperature=.2,kappa=1.,
                               min_gap=1e-9):
    """Explicit spectral Ward partition, shapes 2x3, nondegenerate bands.

    Require T>0; exact or nearly degenerate bands are rejected rather than
    assigning arbitrary Abelian eigenvectors inside a degenerate cluster.
    """
    H=np.asarray(H,complex);V=np.asarray(V,complex);S=np.asarray(S,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError("Hermitian H must be NxN, N>=2")
    n=H.shape[0]
    if V.shape!=(2,n,n) or S.shape!=(3,n,n):
        raise ValueError("V (2,N,N), S (3,N,N) required")
    if temperature<=0 or not np.isfinite(temperature):
        raise ValueError("positive finite temperature required")
    if min_gap<=0 or not np.isfinite(min_gap):
        raise ValueError("positive finite min_gap required")
    for A in (H,*V,*S):
        if not np.all(np.isfinite(A)) or not np.allclose(A,A.conj().T,atol=1e-10,rtol=0):
            raise ValueError("every operator must be finite Hermitian")
    e,U=np.linalg.eigh(H)
    if np.min(np.diff(e))<=min_gap*max(1.,float(np.max(abs(e)))):
        raise ValueError("single-band Ward split requires nondegenerate bands")
    f=_fermi(e,mu,temperature)
    fp=-f*(1.-f)/temperature
    v=np.einsum("in,aij,jm->anm",U.conj(),V,U,optimize=True)
    s=np.einsum("in,bij,jm->bnm",U.conj(),S,U,optimize=True)
    inter=np.zeros((2,3),complex)
    intra=np.zeros((2,3),complex)
    for nn in range(n):
        intra-=kappa*fp[nn]*np.outer(v[:,nn,nn],s[:,nn,nn])
        for mm in range(n):
            if nn==mm:continue
            dd=(f[nn]-f[mm])/(e[nn]-e[mm])
            inter-=kappa*dd*np.outer(v[:,nn,mm],s[:,mm,nn])
    combined=inter+intra
    rho=(U*f)@U.conj().T
    magnetization=np.array([np.trace(rho@b) for b in S])
    return dict(interband=inter,occupation_derivative=intra,
                full_static=combined,minus_kappa_spin_k_gradient=combined.copy(),
                spin_magnetization=magnetization,energies=e,
                fermi_occupations=f,fermi_derivative=fp,
                gap=float(np.min(np.diff(e))))


def equilibrium_current(H,V,*,mu=0.,temperature=.2):
    """Physical Bloch current Tr[f(H) (dH/dk)] at one fixed k."""
    H=np.asarray(H,complex);V=np.asarray(V,complex)
    e,U=np.linalg.eigh(H)
    f=_fermi(e,mu,temperature)
    rho=(U*f)@U.conj().T
    return np.array([np.trace(rho@x) for x in V])


def zeeman_finite_difference(H,V,S,*,mu=0.,temperature=.2,kappa=1.,
                             field_step=2e-5):
    """Independent field-difference oracle: d/dB Tr[f(H-kappa B S)V]."""
    if field_step<=0 or not np.isfinite(field_step):
        raise ValueError("positive field_step required")
    out=np.zeros((2,3),complex)
    for b in range(3):
        plus=equilibrium_current(H-kappa*field_step*S[b],V,
                                 mu=mu,temperature=temperature)
        minus=equilibrium_current(H+kappa*field_step*S[b],V,
                                  mu=mu,temperature=temperature)
        out[:,b]=(plus-minus)/(2*field_step)
    return out


def spin_density(H,S,*,mu=0.,temperature=.2):
    e,U=np.linalg.eigh(H)
    rho=(U*_fermi(e,mu,temperature))@U.conj().T
    return np.array([np.trace(rho@s) for s in S])


def spin_k_gradient_fd(model,k,*,mu=0.,temperature=.2,
                       k_step=2e-5,kappa=1.):
    """Independent finite-k derivative -kappa d_k Tr[f(H) S]."""
    k=np.asarray(k,float)
    if k.shape!=(2,) or k_step<=0:raise ValueError("k and step invalid")
    out=np.zeros((2,3),complex)
    for a in range(2):
        step=np.zeros(2);step[a]=k_step
        hp,_,sp=model(k+step)
        hm,_,sm=model(k-step)
        out[a]=-kappa*(spin_density(hp,sp,mu=mu,temperature=temperature)-
                       spin_density(hm,sm,mu=mu,temperature=temperature))/(2*k_step)
    return out


def static_grid_audit(group,*,order=12,method="midpoint",
                      mu=2.3,temperature=.3,kappa=1.):
    """Unprojected static interband/intraband BZ means and periodic-torus guard.

    Returns the exact discrete telescoping derivative separately. It is NOT
    an assertion that a finite-order spectral quadrature has converged.
    """
    from .primitive_bz import primitive_cell_quadrature
    from .hexagonal_bz import hexagonal_symmetrized_model
    from .group_microscopic import build_group
    from .hexagonal_bz import INV_BASIS
    if method!="midpoint":
        raise ValueError("periodic telescoping gate uses midpoint mesh")
    ks,weights,area=primitive_cell_quadrature(order,method=method)
    ops=build_group(group)
    inter=np.zeros((2,3),complex)
    intra=np.zeros((2,3),complex)
    full=np.zeros((2,3),complex)
    mag=np.zeros((order,order,3),complex)
    for idx,(k,w) in enumerate(zip(ks,weights)):
        h,v,s=hexagonal_symmetrized_model(k,ops)
        r=static_zeeman_current_ward(h,v,s,mu=mu,
                         temperature=temperature,kappa=kappa)
        inter+=w*r["interband"]
        intra+=w*r["occupation_derivative"]
        full+=w*r["full_static"]
        mag[np.unravel_index(idx,(order,order))]=r["spin_magnetization"]
    # Exact periodic nearest-neighbor discrete divergence on fractional torus.
    dq=2*np.pi/order
    gradq=np.stack([(np.roll(mag,-1,axis=a)-np.roll(mag,1,axis=a))/(2*dq)
                   for a in range(2)])
    gradcart=np.einsum("ia,ibcm->abcm",INV_BASIS,gradq,optimize=True)
    telescope=-kappa*np.mean(gradcart,axis=(1,2))
    return dict(interband=inter,occupation_derivative=intra,
                full_static=full,telescoping_derivative=telescope,
                order=order,nodes=len(ks),bz_area=area,
                status="periodic no-current identity; BZ quadrature error NOT bounded")
