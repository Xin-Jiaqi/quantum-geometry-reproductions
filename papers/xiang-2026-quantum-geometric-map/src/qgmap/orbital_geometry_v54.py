"""v54 independently verified multiband Bloch *geometry*, NOT transport.

For nondegenerate fixed, orthonormal, k-independent orbital basis:
r_nm^i=i (U† d_iH U)_nm/(E_m-E_n), q_nm^ij=r_nm^i r_mn^j,
g_nm^ij=Re q_nm^ij, Omega_nm^xy=-2 Im q_nm^xy,
m_tilde_n^xy=1/2 sum_m (E_n-E_m) Omega_nm^xy.
No charges, hbar, magnetic field, or SI units are inserted.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class Geometry:
    energies: np.ndarray
    metric_pairs: np.ndarray
    curvature_pairs: np.ndarray
    metric: np.ndarray
    curvature: np.ndarray
    moment_energy: np.ndarray
    moment_pairs: np.ndarray
    min_direct_gap: float
    acceptance: str = 'VERIFIED_NONDEGENERATE_FIXED_ORBITAL_BASIS'

def band_geometry(hamiltonian,velocities,*,gap_tol=1e-8,eigenvectors=None):
    H=np.asarray(hamiltonian,complex)
    V=np.asarray(velocities,complex)
    if H.ndim<2 or H.shape[-1]!=H.shape[-2] or H.shape[-1]<2:
        raise ValueError('Hamiltonian must be (...,nb,nb) with nb>=2')
    nb=H.shape[-1]
    if V.shape!=(2,)+H.shape:
        raise ValueError('velocities shape must be (2,...,nb,nb) in k_x/k_y order')
    if not np.isfinite(H).all() or not np.isfinite(V).all() or not np.isfinite(gap_tol) or gap_tol<=0:
        raise ValueError('finite operators and positive gap_tol required')
    if not np.allclose(H,np.swapaxes(H.conj(),-1,-2),rtol=0,atol=2e-12):
        raise ValueError('H must be Hermitian')
    if not np.allclose(V,np.swapaxes(V.conj(),-1,-2),rtol=0,atol=2e-12):
        raise ValueError('velocity matrices must be Hermitian')
    E,U=np.linalg.eigh(H)
    gap=np.min(np.diff(E,axis=-1))
    if gap<=gap_tol:
        raise ValueError('unresolved band degeneracy; use a non-Abelian projector solver')
    if eigenvectors is not None:
        U=np.asarray(eigenvectors,complex)
        if U.shape!=H.shape or not np.isfinite(U).all():
            raise ValueError('invalid eigenvectors shape/entries')
        eye=np.eye(nb,dtype=complex)
        if not np.allclose(np.swapaxes(U.conj(),-1,-2)@U,eye,rtol=0,atol=2e-9):
            raise ValueError('provided eigenvectors are not orthonormal')
        if not np.allclose(H@U-U*E[...,None,:],0,rtol=0,atol=2e-8):
            raise ValueError('provided eigenvectors are not matched to sorted H eigenvalues')
    Ud=np.swapaxes(U.conj(),-1,-2)
    vb=np.stack([Ud@V[axis]@U for axis in range(2)],axis=-1)
    delta=E[...,None,:]-E[...,:,None]
    diagonal=np.eye(nb,dtype=bool)
    A=1j*vb/np.where(diagonal,1,delta)[...,None]
    A=np.where(diagonal[...,None],0,A)
    q=A[..., :, :, :, None]*np.swapaxes(A,-3,-2)[..., :, :, None, :]
    g=np.real(q)
    om=-2*np.imag(q[...,0,1])
    weighted=.5*(E[...,:,None]-E[...,None,:])*om
    return Geometry(E,g,om,np.sum(g,axis=-3),np.sum(om,axis=-1),
                    np.sum(weighted,axis=-1),weighted,float(gap))

def curvature_moment_product(geom:Geometry):
    om=geom.curvature_pairs
    mom=geom.moment_pairs
    diag=np.sum(om*mom,axis=-1)
    full=geom.curvature*geom.moment_energy
    return {'full':full,'diagonal_m_eq_l':diag,'cross_m_ne_l':full-diag,
            'acceptance':'VERIFIED_PRODUCT_DECOMPOSITION_NOT_COMPLETE_RESPONSE'}

def fermi_derivative(energy,mu,temperature):
    if not np.isfinite(mu) or not np.isfinite(temperature) or temperature<=0:
        raise ValueError('temperature>0 required')
    f=1/(1+np.exp(np.clip((np.asarray(energy)-mu)/temperature,-450,450)))
    return -f*(1-f)/temperature

def moment_curvature_kernel(geom:Geometry,mu,temperature):
    """ONLY -BZmean(sum f'_n Omega_n m_tilde_n). Not the full orbital MNHE."""
    d=curvature_moment_product(geom)
    fp=fermi_derivative(geom.energies,mu,temperature)
    return {name:float(np.mean(np.sum(-fp*arr,axis=-1)))
            for name,arr in d.items() if name in ('full','diagonal_m_eq_l','cross_m_ne_l')}

def periodic_three_band_model(N=24):
    """Periodic 3-orbital H(k), plus analytically matched dH/dkx,dH/dky."""
    if not isinstance(N,int) or N<6:raise ValueError('N>=6 required')
    k=-np.pi+(np.arange(N)+.33)*2*np.pi/N
    x,y=np.meshgrid(k,k,indexing='ij')
    H=np.zeros((N,N,3,3),complex)
    V=np.zeros((2,N,N,3,3),complex)
    H[...,0,0]=-1.3+.23*np.cos(x)
    H[...,1,1]=.3+.3*np.cos(y)
    H[...,2,2]=2.1+.2*np.cos(x+y)
    V[0,...,0,0]=-.23*np.sin(x)
    V[1,...,1,1]=-.3*np.sin(y)
    V[:,...,2,2]=-.2*np.sin(x+y)
    hops=[
       ((0,1),.43*(np.sin(x)+1j*np.sin(y))+.13j*np.cos(x),
        .43*np.cos(x)-.13j*np.sin(x),1j*.43*np.cos(y)),
       ((0,2),.38*(np.sin(y)+1j*np.sin(x)),
        1j*.38*np.cos(x),.38*np.cos(y)),
       ((1,2),.49*(np.sin(x+y)+1j*np.cos(x)),
        .49*(np.cos(x+y)-1j*np.sin(x)),.49*np.cos(x+y))
    ]
    for (i,j),h,vx,vy in hops:
        H[...,i,j]=h;H[...,j,i]=h.conj()
        V[0,...,i,j]=vx;V[0,...,j,i]=vx.conj()
        V[1,...,i,j]=vy;V[1,...,j,i]=vy.conj()
    return k,H,V

def finite_difference_projector_reference(model,kx,ky,step=1e-4):
    """Independent P=[u><u] finite differences: Omega=iTr P[Px,Py].
    Metric via perpendicular du; moment via Im <du_x|(H-E)|du_y>.
    Input model(x,y) must use a fixed orthonormal orbital basis.
    """
    if step<=0 or not np.isfinite(step):raise ValueError('positive step')
    h0=np.asarray(model(kx,ky),complex)
    E,U=np.linalg.eigh(h0)
    def get(x,y):return np.linalg.eigh(model(x,y))[1]
    ux_p=get(kx+step,ky);ux_m=get(kx-step,ky)
    uy_p=get(kx,ky+step);uy_m=get(kx,ky-step)
    omegas=[];moms=[];metric=[]
    for n in range(len(E)):
        u=U[:,n]
        P=np.outer(u,u.conj())
        px=(np.outer(ux_p[:,n],ux_p[:,n].conj())-np.outer(ux_m[:,n],ux_m[:,n].conj()))/(2*step)
        py=(np.outer(uy_p[:,n],uy_p[:,n].conj())-np.outer(uy_m[:,n],uy_m[:,n].conj()))/(2*step)
        omegas.append(float(np.real(1j*np.trace(P@(px@py-py@px)))))
        def align(v):return v*np.exp(-1j*np.angle(np.vdot(u,v)))
        dx=(align(ux_p[:,n])-align(ux_m[:,n]))/(2*step)
        dy=(align(uy_p[:,n])-align(uy_m[:,n]))/(2*step)
        perp=np.eye(len(E))-P
        dx=perp@dx;dy=perp@dy
        metric.append(float(np.real(np.vdot(dx,dy))))
        moms.append(float(np.imag(np.vdot(dx,(h0-E[n]*np.eye(len(E)))@dy))))
    return {'curvature':np.asarray(omegas),'metric_xy':np.asarray(metric),
            'moment_energy':np.asarray(moms)}
