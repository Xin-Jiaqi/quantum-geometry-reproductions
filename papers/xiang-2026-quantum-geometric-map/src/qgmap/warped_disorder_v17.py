"""v17: 2D mirror-compatible periodic lattice completion with scalar white-noise SCBA and ladder.

This independently chosen tight-binding regulator is NOT the PRB's continuum TI surface,
nor a lattice fitting of the paper's static PHE. It necessarily contains extra Dirac cones
and its uniform onsite Zeeman cannot equal a globally periodic k shift.

H0(k) = beta sin kx I + 2alpha(2-cos kx-cos ky)I
      + sin kx sigma_y - sin ky sigma_x
      + lam sin ky (sin^2 ky - 3 sin^2 kx) sigma_z.
Zeeman branch: H=H0-b sigma_x. Shift control: H=H0(kx,ky+b).

The observable is Re[Tr(v_x GR Gamma_y GA+v_y GR Gamma_x GA)/2]/N^2
at fixed energy mu with finite eta, scalar Gaussian onsite disorder within SCBA
and a noncrossing RA current ladder. This is not full Kubo-Streda, microscopic
finite frequency, or the paper's two-color Liouville/RTA coefficient.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .microscopic_ward_v11 import I,SX,SY,SZ,scba,ladder_vertex,greens

@dataclass(frozen=True)
class Lattice:
    n:int=101
    beta:float=.125
    alpha:float=.16
    lam:float=.6
    chemical_potential:float=.55
    eta:float=.08
    disorder_variance:float=.045
    offset_x:float=.31
    offset_y:float=.37
    def __post_init__(self):
        if self.n<12 or self.eta<=0 or self.disorder_variance<0 or abs(self.beta)>=1:
            raise ValueError('Invalid lattice, regulator, disorder or tilt')


def k_grid(p:Lattice):
    kx=-np.pi+2*np.pi*(np.arange(p.n)+p.offset_x)/p.n
    ky=-np.pi+2*np.pi*(np.arange(p.n)+p.offset_y)/p.n
    return np.stack(np.meshgrid(kx,ky,indexing='ij'),axis=-1)


def h_and_vel(k,p:Lattice,b=0.,protocol='onsite'):
    k=np.asarray(k,dtype=float)
    x=k[...,0]; y=k[...,1]
    if protocol=='shift':
        y=y+b
        field=0.
    elif protocol=='onsite':
        field=b
    else:
        raise ValueError('protocol should be onsite or shift')
    sx,sy=np.sin(x),np.sin(y)
    cx,cy=np.cos(x),np.cos(y)
    scalar=p.beta*sx+2*p.alpha*(2-cx-cy)
    W=sy*(sy*sy-3*sx*sx)
    H=scalar[...,None,None]*I + sx[...,None,None]*SY -sy[...,None,None]*SX + (p.lam*W)[...,None,None]*SZ - field*SX
    Wx=-6*sy*sx*cx
    Wy=3*cy*(sy*sy-sx*sx)
    vx=(p.beta*cx+2*p.alpha*sx)[...,None,None]*I+cx[...,None,None]*SY+(p.lam*Wx)[...,None,None]*SZ
    vy=(2*p.alpha*sy)[...,None,None]*I-cy[...,None,None]*SX+(p.lam*Wy)[...,None,None]*SZ
    return H,vx,vy


def matrix_components(k,p:Lattice):
    """Exact mirror M_y: U=i sigma_y, (kx,ky,b)->(kx,-ky,-b)."""
    return h_and_vel(k,p)


def ra_observable(p:Lattice,b=0.,protocol='onsite',vertex=True):
    k=k_grid(p)
    h,vx,vy=h_and_vel(k,p,b,protocol)
    sig,gr,ga,niter=scba(h,p)
    if vertex:
        gy,Cy,resy=ladder_vertex(gr,ga,vy,p.disorder_variance)
        gx,Cx,resx=ladder_vertex(gr,ga,vx,p.disorder_variance)
    else:
        gx,gy=vx,vy
        Cx=Cy=np.zeros((2,2),complex)
        resx=resy=0.
    xy=np.mean(np.trace(vx@gr@gy@ga,axis1=-2,axis2=-1))
    yx=np.mean(np.trace(vy@gr@gx@ga,axis1=-2,axis2=-1))
    return dict(symmetric=float(np.real((xy+yx)/2)),xy=float(np.real(xy)),yx=float(np.real(yx)),
                vertex_residual=float(max(resx,resy)),scba_iterations=niter,
                self_energy=sig,vertex_x=Cx,vertex_y=Cy,
                scba_residual=float(np.max(abs(sig-p.disorder_variance*gr.mean((0,1))))))


def slope(p:Lattice,step=.025,protocol='onsite',vertex=True):
    a=ra_observable(p,step,protocol,vertex)
    b=ra_observable(p,-step,protocol,vertex)
    return dict(value=(a['symmetric']-b['symmetric'])/(2*step),
                plus=a['symmetric'],minus=b['symmetric'],
                self_energy_change=float(np.max(abs(a['self_energy']-b['self_energy']))),
                vertex_residual=max(a['vertex_residual'],b['vertex_residual']),
                scba_residual=max(a['scba_residual'],b['scba_residual']))


def density_ward(p:Lattice,b=0.,protocol='onsite'):
    h,_,_=h_and_vel(k_grid(p),p,b,protocol)
    sig,gr,ga,niter=scba(h,p)
    vertex,constant,res=ladder_vertex(gr,ga,np.broadcast_to(I,h.shape),p.disorder_variance)
    from_sigma=I-(sig-sig.conj().T)/(2j*p.eta)
    return dict(abs_error=float(np.max(abs(I+constant-from_sigma)),
                residual=res,iterations=niter,
                causality_max_eigenvalue=float(np.max(np.linalg.eigvalsh((sig-sig.conj().T)/(2j)))))
