"""Round27: exact finite-q longitudinal U(1) Ward on a finite 2D SOC torus.

Finite Hilbert space: one spin-1/2 orbital on each site, hbar=1, charge
absorbed into the onsite scalar potential and Peierls phase. All energy and
frequency arguments are in the same natural energy unit. A finite torus is
NOT a bulk metallic KME or a conserving collision/bath model.

Time convention e^{-i z t}, z=omega+i eta with eta>0. For a real diagonal
site-phase generator L, a pure gauge U=exp(i lambda(t)L) gives
H_g = U H U^dagger - dot(lambda)L. First-order vector drive K=i[L,H],
scalar drive V_s=i z L; V_g=K+V_s=i(z-ad_H)L.

rho_g response = i[L,rho_0]; every local density unchanged; physical bond
current response includes i[L,J] contact and vanishes for pure gauge.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

_X=np.array([[0,1],[1,0]],complex)
_Y=np.array([[0,-1j],[1j,0]],complex)
_Z=np.diag([1.,-1.]).astype(complex)
_I=np.eye(2,dtype=complex)

@dataclass
class FiniteTorus:
    H: np.ndarray
    bonds: tuple
    projectors: np.ndarray
    shape: tuple

    def check(self):
        h=np.asarray(self.H,complex)
        if h.ndim!=2 or h.shape[0]!=h.shape[1] or h.shape[0]%2 or not np.all(np.isfinite(h)):
            raise ValueError('finite even-dimensional Hamiltonian required')
        if not np.allclose(h,h.conj().T,atol=1e-12,rtol=0):
            raise ValueError('Hermitian H required')
        if len(self.shape)!=2 or min(self.shape)<3 or h.shape[0]!=2*np.prod(self.shape):
            raise ValueError('finite 2D torus with at least 3 sites along each axis')
        if self.projectors.shape!=(np.prod(self.shape),*h.shape):
            raise ValueError('site projectors required')
        if len(self.bonds)!=2*np.prod(self.shape):
            raise ValueError('positive x/y directed bonds required')
        return self


def soc_torus(nx=3,ny=3,*,t_x=.83,t_y=.66,alpha=.23,beta=.17,m=.37):
    """Explicit spin-orbit lattice H: x: t_x I+i alpha sigma_y;
    y: t_y I-i beta sigma_x; onsite m sigma_z + weak scalar modulation.
    """
    if not all(isinstance(v,int) and v>=3 for v in (nx,ny)):
        raise ValueError('minimum 3 by 3 torus')
    if not all(np.isfinite(x) for x in (t_x,t_y,alpha,beta,m)):
        raise ValueError('finite hopping/SOC/Zeeman model required')
    sites=nx*ny;n=2*sites
    H=np.zeros((n,n),complex);P=np.zeros((sites,n,n),complex);bonds=[]
    Tx=t_x*_I+1j*alpha*_Y
    Ty=t_y*_I-1j*beta*_X
    ix=lambda x,y:(y%ny)*nx+x%nx
    for y in range(ny):
        for x in range(nx):
            i=ix(x,y);sl=slice(2*i,2*i+2)
            H[sl,sl]+=m*_Z+(.07*np.cos(2*np.pi*x/nx)+.04*np.sin(2*np.pi*y/ny))*_I
            P[i,sl,sl]=_I
            for axis,j,T in [('x',ix(x+1,y),Tx),('y',ix(x,y+1),Ty)]:
                sj=slice(2*j,2*j+2)
                H[sl,sj]+=T;H[sj,sl]+=T.conj().T
                bonds.append((i,j,axis))
    return FiniteTorus(H,tuple(bonds),P,(nx,ny)).check()


def site_phase(torus,*,qx=1,qy=1,offset=.21):
    """Real single nonzero lattice Fourier mode (cosine at q and -q).

    On the Nx by Ny torus, q=(2pi qx/Nx,2pi qy/Ny). The site phase is
    defined in real space, no small-q/continuum gradient approximation.
    """
    torus.check();nx,ny=torus.shape
    if not (isinstance(qx,int) and isinstance(qy,int) and
            0<=qx<nx and 0<=qy<ny and (qx,qy)!=(0,0)):
        raise ValueError('nonzero finite lattice Fourier q required')
    if not np.isfinite(offset):raise ValueError('finite gauge phase required')
    vals=np.array([np.cos(2*np.pi*(qx*x/nx+qy*y/ny)+offset)
                   for y in range(ny) for x in range(nx)])
    return np.diag(np.repeat(vals,2)).astype(complex)


def bond_currents(torus,H=None):
    """J_{i->j}=i(H_ij |i><j| - h.c.), outward from i.
    The same definition must be used for field-dependent Peierls H.
    """
    torus.check();h=torus.H if H is None else np.asarray(H,complex)
    if h.shape!=torus.H.shape or not np.allclose(h,h.conj().T,rtol=0,atol=1e-10):
        raise ValueError('matched Hermitian field Hamiltonian required')
    out=[]
    for i,j,_ in torus.bonds:
        z=np.zeros_like(h)
        z[2*i:2*i+2,2*j:2*j+2]=1j*h[2*i:2*i+2,2*j:2*j+2]
        z[2*j:2*j+2,2*i:2*i+2]=-1j*h[2*j:2*j+2,2*i:2*i+2]
        out.append(z)
    return np.asarray(out)


def divergence(torus,J):
    if np.asarray(J).shape!=(len(torus.bonds),*torus.H.shape):
        raise ValueError('one matrix per directed positive bond required')
    out=np.zeros((len(torus.projectors),*torus.H.shape),complex)
    for a,(i,j,_) in enumerate(torus.bonds):
        out[i]+=J[a]
        out[j]-=J[a]
    return out


def fermi_density(H,mu=.2,temperature=.3):
    if not np.isfinite(mu) or not np.isfinite(temperature) or temperature<=0:
        raise ValueError('positive finite T required; collisionless grand canonical')
    e,u=np.linalg.eigh(H)
    f=1/(1+np.exp(np.clip((e-mu)/temperature,-700,700)))
    return (u*f)@u.conj().T


def linear_spectral(H,rho,V,z):
    """Independent eigenbasis Kubo density response (z-ad_H)^-1[V,rho]."""
    H=np.asarray(H,complex);V=np.asarray(V,complex);rho=np.asarray(rho,complex)
    z=complex(z)
    if z.imag<=0 or not np.isfinite(z):
        raise ValueError('retarded Im z >0 required')
    if V.shape!=H.shape or rho.shape!=H.shape:
        raise ValueError('matched H,V,rho required')
    e,u=np.linalg.eigh(H)
    vb=u.conj().T@V@u
    rb=u.conj().T@rho@u
    comm=vb@rb-rb@vb
    diff=e[:,None]-e[None,:]
    return u@(comm/(z-diff))@u.conj().T


def linear_liouville(H,rho,V,z):
    """Independent full orbital-basis N^2 Liouville inversion; small torus only."""
    H=np.asarray(H,complex);n=len(H);z=complex(z)
    if z.imag<=0:raise ValueError('retarded Im z >0 required')
    ad=np.kron(np.eye(n),H)-np.kron(H.T,np.eye(n))
    C=V@rho-rho@V
    out=np.linalg.solve(z*np.eye(n*n)-ad,C.reshape(-1,order='F'))
    return out.reshape((n,n),order='F')


def pure_gauge_response(torus,L,*,z=.42+.09j,mu=.2,temperature=.3):
    """Split scalar/vector response and return full current contact and densities."""
    torus.check();h=torus.H;L=np.asarray(L,complex);z=complex(z)
    if L.shape!=h.shape or not np.allclose(L,L.conj().T,atol=1e-12) or not np.allclose(L,np.diag(np.diag(L)),atol=1e-12):
        raise ValueError('real diagonal position gauge required')
    if z.imag<=0:raise ValueError('retarded Im z >0 required')
    rho=fermi_density(h,mu,temperature)
    K=1j*(L@h-h@L)
    Vscalar=1j*z*L
    dr_vector=linear_spectral(h,rho,K,z)
    dr_scalar=linear_spectral(h,rho,Vscalar,z)
    dr_gauge=dr_vector+dr_scalar
    dr_exact=1j*(L@rho-rho@L)
    J=bond_currents(torus)
    J_contact=1j*(L[None]@J-J@L[None])
    density=lambda dr:np.einsum('aij,ji->a',torus.projectors,dr,optimize=True)
    j_state=lambda dr:np.einsum('aij,ji->a',J,dr,optimize=True)
    j_contact=np.einsum('aij,ji->a',J_contact,rho,optimize=True)
    return dict(rho=rho,L=L,z=z,K=K,Vscalar=Vscalar,
           dr_vector=dr_vector,dr_scalar=dr_scalar,dr_gauge=dr_gauge,
           dr_exact=dr_exact,density_vector=density(dr_vector),
           density_scalar=density(dr_scalar),density_gauge=density(dr_gauge),
           current_state_vector=j_state(dr_vector),
           current_state_scalar=j_state(dr_scalar),
           current_contact=j_contact,
           current_gauge=j_state(dr_gauge)+j_contact,
           bond_J=J,bond_J_contact=J_contact)


def continuity_residual(torus,dr,V,rho,z,J_contact=None):
    """Finite-q continuity: div(delta J) = i z delta n.
    Input J_contact is explicit first-order field derivative of bond current.
    """
    js=bond_currents(torus)
    delta=np.einsum('bij,ji->b',js,dr,optimize=True)
    if J_contact is not None:delta+=np.einsum('bij,ji->b',J_contact,rho,optimize=True)
    div=np.zeros(len(torus.projectors),complex)
    for b,(i,j,_) in enumerate(torus.bonds):div[i]+=delta[b];div[j]-=delta[b]
    dn=np.einsum('aij,ji->a',torus.projectors,dr,optimize=True)
    return div-1j*z*dn


def finite_amplitude_gauge(torus,L,amplitude):
    """Independent finite Peierls phase via exact unitary, no linear truncation."""
    torus.check();L=np.asarray(L,complex)
    if not np.allclose(L,np.diag(np.diag(L)),rtol=0,atol=1e-12):
        raise ValueError('diagonal phase needed')
    U=np.diag(np.exp(1j*amplitude*np.diag(L)))
    return U,U@torus.H@U.conj().T
