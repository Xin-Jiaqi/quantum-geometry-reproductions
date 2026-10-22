"""First independent checkpoint for Cui et al. (2026) Eqs. (2)-(15).

Units: dimensionless k, energy; H=d(k) dot sigma. The Pauli operator is
specified as sigma (not hbar sigma/2). Outputs are geometric, NOT SI transport.
No nondegenerate band formula is evaluated at |d|=0.
"""
from __future__ import annotations
import numpy as np

PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)

def model(kx:float,ky:float,m:float=0.,vx:float=1.,vy:float=1.):
    d=np.array([vx*kx,vy*ky,m],float)
    deriv=np.array([[vx,0.,0.],[0.,vy,0.],[0.,0.,0.]],float)
    return d,deriv

def _cross_geom(T):
    """Eq. (3), using a,b ordered x,y,z; Omega_c=(1/2)eps_cab Omega_ab."""
    sym=(T+T.T)/2
    anti=T-T.T
    gN=sym.real; gA=sym.imag
    # Omega^{ab} = -Re/Im[T^{ab}-T^{ba}]; vectors Omega_x=Omega^{yz}, etc.
    AN=-anti.imag; AA=-anti.real
    vN=np.array([AN[1,2],AN[2,0],AN[0,1]])
    vA=np.array([AA[1,2],AA[2,0],AA[0,1]])
    return dict(gN=gN,gA=gA,omegaN=vN,omegaA=vA)

def qgt_eigen(d,dk,band:int=-1,spin_basis=None):
    """Original interband Eq.(2) using eigensystem and r_nm^a.

    Optional spin_basis rotates physical Pauli operators together with H and dH,
    ensuring invariance under basis changes (not phase only).
    """
    E=float(np.linalg.norm(d))
    if E<=1e-12: raise ValueError('Degenerate Dirac point; interband formula invalid')
    if band not in (-1,1): raise ValueError('band must be -1 or 1')
    H=np.einsum('a,aij->ij',d,PAULI)
    dH=np.einsum('xa,aij->xij',dk,PAULI)
    sig=PAULI.copy()
    if spin_basis is not None:
        U=np.asarray(spin_basis,dtype=complex)
        H=U@H@U.conj().T
        dH=np.einsum('ij,ajk,kl->ail',U,dH,U.conj().T)
        sig=np.einsum('ij,ajk,kl->ail',U,sig,U.conj().T)
    es,us=np.linalg.eigh(H)
    n=0 if band==-1 else 1; j=1-n
    en,em=es[n],es[j]
    ket=us[:,n]; other=us[:,j]
    r=np.array([1j*np.vdot(ket,v@other)/(em-en) for v in dH])
    s=np.array([np.vdot(other,o@ket) for o in sig])
    T=np.outer(r,s)
    return dict(T=T,**_cross_geom(T))

def qgt_projector(d,dk,band:int=-1):
    """Independent projector/trace Eq.(36) with spectral projectors (no gauge)."""
    E=float(np.linalg.norm(d))
    if E<=1e-12: raise ValueError('Degenerate Dirac point')
    if band not in (-1,1): raise ValueError('band must be -1 or 1')
    H=np.einsum('a,aij->ij',d,PAULI)
    Pn=(np.eye(2)+band*H/E)/2
    Pm=np.eye(2)-Pn
    dH=np.einsum('xa,aij->xij',dk,PAULI)
    denom=-2*band*E
    T=np.array([[1j*np.trace(Pn@dH[a]@Pm@PAULI[b])/denom for b in range(3)] for a in range(3)])
    return dict(T=T,**_cross_geom(T))

def qgt_analytic(d,dk,band:int=-1):
    """Independent compact Eq.(53) and Eq.(3)."""
    E=float(np.linalg.norm(d))
    if E<=1e-12: raise ValueError('Degenerate Dirac point')
    n=d/E
    cross=np.cross(n,dk)
    perpendicular=dk-np.einsum('a,xa->x',n,dk)[:,None]*n
    T=cross/(2*E)-1j*perpendicular/(2*band*E)
    return dict(T=T,**_cross_geom(T))

def omega_analytic_massive(kx,ky,m=0.):
    """Main Eq.(6)-(8) for d=(kx,ky,m), returned [Omega_x,Omega_y,Omega_z]."""
    D=kx*kx+ky*ky+m*m
    if D<=0: raise ValueError('Singular Dirac point')
    return -np.array([kx,ky,2*m],float)/(2*D)

def winding_field(kx,ky,vx=1.,vy=1.):
    """Analytic phase-gradient for planar d=(vx*kx,vy*ky,0)."""
    D=(vx*kx)**2+(vy*ky)**2
    if D<=0: raise ValueError('Singular Dirac point')
    return np.array([-vx*vy*ky,vx*vy*kx])/D

def flux_circle(radius,m=0.,points=512):
    """Numerical normal line integral of Eq.(8) on a circular contour."""
    th=2*np.pi*np.arange(points)/points
    xy=np.array([radius*np.cos(th),radius*np.sin(th)]).T
    ds=np.roll(xy,-1,axis=0)-xy
    mids=(xy+np.roll(xy,-1,axis=0))/2
    dr=np.array([ds[:,1],-ds[:,0]]).T
    F=np.array([omega_analytic_massive(*xy,m)[:2] for xy in mids])
    return float(np.sum(F*dr))

def flux_ellipse(a,b,vx=1.,vy=1.,points=4096):
    """Flux of planar anisotropic Dirac, normal dl=(dy,-dx), CCW."""
    th=2*np.pi*np.arange(points)/points
    xy=np.array([a*np.cos(th),b*np.sin(th)]).T
    dxy=np.roll(xy,-1,axis=0)-xy
    mids=(xy+np.roll(xy,-1,axis=0))/2
    rotated_normal=np.array([dxy[:,1],-dxy[:,0]]).T
    w=np.array([winding_field(*x,vx=vx,vy=vy) for x in mids])
    F=np.array([-w[:,1],w[:,0]]).T/2
    return float(np.sum(F*rotated_normal))

def winding_ellipse(a,b,vx=1.,vy=1.,points=4096):
    th=2*np.pi*np.arange(points)/points
    xy=np.array([a*np.cos(th),b*np.sin(th)]).T
    dm=np.roll(xy,-1,axis=0)-xy
    mid=(xy+np.roll(xy,-1,axis=0))/2
    w=np.array([winding_field(*x,vx=vx,vy=vy) for x in mid])
    return float(np.sum(w*dm)/(2*np.pi))

def berry_phase_wilson_circle(radius,m=0.,points=512,band=-1):
    """Independent Wilson-loop Berry phase with spinor overlaps (mod 2pi)."""
    th=2*np.pi*np.arange(points)/points
    us=[]
    for t in th:
        d,_=model(radius*np.cos(t),radius*np.sin(t),m)
        _,u=np.linalg.eigh(np.einsum('a,aij->ij',d,PAULI))
        us.append(u[:,0 if band==-1 else 1])
    return float(-np.angle(np.prod([np.vdot(us[i],us[(i+1)%points])/abs(np.vdot(us[i],us[(i+1)%points])) for i in range(points)])))

def flux_circle_from_qgt(radius,m=0.,points=512,band=-1):
    """Integrate the raw Zeeman QGT eigenvector matrix elements, not Eq.(8)."""
    th=2*np.pi*np.arange(points)/points
    xy=np.array([radius*np.cos(th),radius*np.sin(th)]).T
    ds=np.roll(xy,-1,axis=0)-xy
    mids=(xy+np.roll(xy,-1,axis=0))/2
    rotated_normals=np.array([ds[:,1],-ds[:,0]]).T
    omegas=[]
    for kx,ky in mids:
        d,dk=model(kx,ky,m)
        omegas.append(qgt_eigen(d,dk,band)['omegaA'][:2])
    return float(np.sum(np.asarray(omegas)*rotated_normals))

def divergence_disk_massive(radius,m,points=10000):
    """2D surface integral of div OmegaA = -m²/(r²+m²)², independent Gauss route."""
    if abs(m)<=1e-14: raise ValueError('massless divergence is a distribution, not an ordinary function')
    r=(np.arange(points)+.5)*radius/points
    return float(np.sum(-m*m/(r*r+m*m)**2*2*np.pi*r*radius/points))
