"""A disorder-self-energy / ladder-vertex benchmark for a magnetic momentum-shift Ward identity.

This is a deliberately INDEPENDENT model audit, NOT a claim to reproduce the
Xiang et al. relaxation-time coefficient. At fixed Fermi energy it computes the
RA Fermi-surface bubble (not a full Kubo-Streda conductivity), with a scalar
momentum-independent impurity variance and a retarded self-consistent Born
self-energy plus a noncrossing ladder current vertex.

Two field protocols:
    shift:   H(k,B)=H0(kx,ky+B), J_a=partial_k_a H(k,B)
    zeeman:  H(k,B)=H0(k) - B sigma_x, J_a=partial_k_a H0(k)

The first model has exact k-space shift covariance on the full periodic BZ:
Sigma(z,B)=Sigma(z,0), Lambda_a(k,B)=Lambda_a(kx,ky+B,0), so the
integrated RA current response is B independent for a translation-invariant
regulator. The second field is an ONSITE Zeeman perturbation; it breaks this
identity for the Wilson lattice completion. Neither is asserted to describe
microscopic scattering in the published continuum surface model.

The nonzero external eta is a retarded resolvent regulator. Finite eta means
this is *not* a proof of the fully charge-conserving finite-q Ward identity.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
BASIS=np.array([I,SX,SY,SZ])

@dataclass(frozen=True)
class WardParameters:
    n:int=61
    beta:float=.125
    wilson_mass:float=.35
    chemical_potential:float=.22
    eta:float=.08
    disorder_variance:float=.045
    offset_x:float=.37
    offset_y:float=.29

    def __post_init__(self):
        if self.n<8 or self.eta<=0 or self.disorder_variance<0 or not 0<=abs(self.beta)<1:
            raise ValueError('Invalid self-energy / mesh parameters')


def momentum_grid(p:WardParameters):
    kx=-np.pi+2*np.pi*(np.arange(p.n)+p.offset_x)/p.n
    ky=-np.pi+2*np.pi*(np.arange(p.n)+p.offset_y)/p.n
    return np.stack(np.meshgrid(kx,ky,indexing='ij'),axis=-1)


def h0(k,p):
    x,y=k[...,0],k[...,1]
    return (p.beta*np.sin(x)[...,None,None]*I
            +np.sin(x)[...,None,None]*SY
            -np.sin(y)[...,None,None]*SX
            +p.wilson_mass*(2-np.cos(x)-np.cos(y))[...,None,None]*SZ)


def band_model(k,p,B=0.,coupling='shift'):
    if coupling=='shift':
        k2=np.array(k,copy=True)
        k2[...,1]+=B
        H=h0(k2,p)
        ky=k2[...,1]
    elif coupling=='zeeman':
        H=h0(k,p)-B*SX
        ky=k[...,1]
    else:
        raise ValueError('Choose explicit shift or Zeeman perturbation')
    kx=k[...,0]
    vx=(np.cos(kx)[...,None,None]*(p.beta*I+SY)
        +p.wilson_mass*np.sin(kx)[...,None,None]*SZ)
    vy=(-np.cos(ky)[...,None,None]*SX
        +p.wilson_mass*np.sin(ky)[...,None,None]*SZ)
    return H,vx,vy


def inv2(m):
    """Vectorized independent adjugate, no band eigenvectors / gauge choice."""
    a,b=m[...,0,0],m[...,0,1]
    c,d=m[...,1,0],m[...,1,1]
    det=a*d-b*c
    out=np.empty_like(m)
    out[...,0,0]=d/det
    out[...,0,1]=-b/det
    out[...,1,0]=-c/det
    out[...,1,1]=a/det
    return out


def greens(h,z,sigma):
    return inv2(z*I-h-sigma)


def scba(h,p,max_iterations=300,tolerance=2e-12):
    """Sigma^R=W <(mu+i eta-H-Sigma^R)^-1> with stabilizing mixing."""
    sigma=np.zeros((2,2),complex)
    z=p.chemical_potential+1j*p.eta
    for it in range(max_iterations):
        gr=greens(h,z,sigma)
        cand=p.disorder_variance*gr.mean(axis=(0,1))
        new=.6*cand+.4*sigma
        if np.max(np.abs(new-sigma))<tolerance:
            sigma=new
            break
        sigma=new
    else:
        raise RuntimeError('SCBA failed to converge')
    gr=greens(h,z,sigma)
    ga=np.swapaxes(gr.conj(),-1,-2)
    return sigma,gr,ga,it+1


def ladder_vertex(gr,ga,v,impurity_variance):
    """Gamma(k)=v(k)+C, C=W <GR(v(k)+C)GA>. Full 4x4 Pauli basis."""
    W=impurity_variance
    phi=np.mean(gr@v@ga,axis=(0,1))
    mat=np.empty((4,4),complex)
    rhs=np.empty(4,complex)
    for a,t in enumerate(BASIS):
        rhs[a]=.5*np.trace(t@phi)
        for b,s in enumerate(BASIS):
            kernel=np.mean(gr@s@ga,axis=(0,1))
            mat[a,b]=.5*np.trace(t@kernel)
    coeff=np.linalg.solve(np.eye(4)-W*mat,W*rhs)
    C=np.einsum('a,aij->ij',coeff,BASIS)
    residual=C-W*np.mean(gr@(v+C)@ga,axis=(0,1))
    return v+C,C,float(np.max(abs(residual)))


def ra_kernel(k,p,B=0.,coupling='shift',vertex_corrections=True):
    """Symmetric RA Fermi-surface kernel in arbitrary dimensionless units."""
    h,vx,vy=band_model(k,p,B,coupling)
    sigma,gr,ga,niter=scba(h,p)
    if vertex_corrections:
        vertex_y,Cy,ry=ladder_vertex(gr,ga,vy,p.disorder_variance)
        vertex_x,Cx,rx=ladder_vertex(gr,ga,vx,p.disorder_variance)
    else:
        vertex_y,Cy,ry=vy,np.zeros((2,2),complex),0.
        vertex_x,Cx,rx=vx,np.zeros((2,2),complex),0.
    # RA transport bubble, not a full Kubo-Streda conductivity.
    xy=np.mean(np.trace(vx@gr@vertex_y@ga,axis1=-2,axis2=-1))
    yx=np.mean(np.trace(vy@gr@vertex_x@ga,axis1=-2,axis2=-1))
    return dict(symmetric=float(np.real(.5*(xy+yx))),
                xy=float(np.real(xy)),yx=float(np.real(yx)),
                scba_iterations=niter, selfenergy=sigma, vertex_y=Cy,
                vertex_x=Cx,vertex_residual=max(ry,rx),
                spectral_imaginary_part=float(np.max(np.linalg.eigvalsh((sigma-sigma.conj().T)/(2j)))))


def response_slope(k,p,step=.02,coupling='shift',vertex_corrections=True):
    plus=ra_kernel(k,p,+step,coupling,vertex_corrections)
    minus=ra_kernel(k,p,-step,coupling,vertex_corrections)
    zero=ra_kernel(k,p,0.,coupling,vertex_corrections)
    return dict(slope=(plus['symmetric']-minus['symmetric'])/(2*step),
                at_zero=zero['symmetric'],
                at_plus=plus['symmetric'],at_minus=minus['symmetric'],
                selfenergy_B_delta=float(np.max(abs(plus['selfenergy']-minus['selfenergy']))),
                vertex_residual=max(plus['vertex_residual'],minus['vertex_residual']),
                iteration_max=max(plus['scba_iterations'],minus['scba_iterations']))


def density_ward_identity(k,p,B=0.,coupling='shift'):
    """SCBA scalar-disorder q=0 RA density-vertex Ward check (finite eta).

    Gamma_0 = I - (Sigma_R - Sigma_A)/(2i eta),
    while the *independent* Bethe-Salpeter density vertex is
    Gamma_0 = I + W <GR Gamma_0 GA>.
    Returns the max matrix-element discrepancy and separate residuals.
    Finite q, the strict eta->0 limit and all crossed diagrams are NOT tested.
    """
    H,_,_=band_model(k,p,B,coupling)
    sig,gr,ga,its=scba(H,p)
    sr_sa=sig-sig.conj().T
    from_selfenergy=I - sr_sa/(2j*p.eta)
    vertex,correction,residual=ladder_vertex(gr,ga,np.broadcast_to(I,H.shape),p.disorder_variance)
    from_ladder=I+correction
    return dict(max_abs_diff=float(np.max(abs(from_ladder-from_selfenergy))),
                ladder_residual=residual,
                selfenergy_residual=float(np.max(abs(sig-p.disorder_variance*gr.mean((0,1))))),
                gamma_from_ward=from_selfenergy,
                gamma_from_ladder=from_ladder,scba_iterations=its)
