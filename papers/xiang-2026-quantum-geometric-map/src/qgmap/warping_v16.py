"""v16 independently defined TI surface model extension, not the PRB Eq.(10).

H(k,b) = [beta*kx+alpha*(kx**2+ky**2)] I + kx*sigma_y
         -(ky+b)*sigma_x + lam*W(kx,ky)*sigma_z.

Here b is the Zeeman momentum scale, not tesla, and e=hbar=v_F=tau=1.
All responses are 2D static re-thermalized constant-tau Drude, not the
published two-regulator Liouville/RTA prescription, and use model units.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np

@dataclass(frozen=True)
class Surface:
    beta: float=.125
    alpha: float=0.
    lam: float=0.
    warp_axis: str="mirror_y"
    mu: float=.55
    temp: float=.065
    radius: float=3.2

    def __post_init__(self):
        if self.temp <= 0 or self.radius <= 0 or abs(self.beta)>=1 or self.warp_axis not in ("mirror_y","mirror_x"):
            raise ValueError('Requires temp,radius>0 and |tilt beta|<1')

@lru_cache(maxsize=30)
def grid(nrad=140,nang=256,radius=3.2):
    from numpy.polynomial.legendre import leggauss
    if nrad < 2 or nang < 8 or radius <= 0:
        raise ValueError('grid resolution or cutoff invalid')
    x,w=leggauss(nrad)
    r=radius*(x+1)/2
    wr=w*radius/2*r/(2*np.pi*nang)
    theta=np.arange(nang)*2*np.pi/nang
    kx=(r[:,None]*np.cos(theta)[None,:])
    ky=(r[:,None]*np.sin(theta)[None,:])
    return kx,ky,np.broadcast_to(wr[:,None],kx.shape)


def _elements(p:Surface,kx,ky,b=0.,shifted=False):
    if shifted:
        yphys=ky-b
        dx=-ky
    else:
        yphys=ky
        dx=-(ky+b)
    dy=kx
    if p.warp_axis=="mirror_y":
        dz=p.lam*(yphys**3-3*yphys*kx**2)
        dz_x=-6*p.lam*kx*yphys
        dz_y=3*p.lam*(yphys*yphys-kx*kx)
        dzb=-dz_y
        dzxb=6*p.lam*kx
        dzyb=-6*p.lam*yphys
    else:
        dz=p.lam*(kx**3-3*kx*yphys**2)
        dz_x=3*p.lam*(kx*kx-yphys*yphys)
        dz_y=-6*p.lam*kx*yphys
        dzb=-dz_y
        dzxb=6*p.lam*yphys
        dzyb=6*p.lam*kx
    c=p.beta*kx+p.alpha*(kx*kx+yphys*yphys)
    c_x=p.beta+2*p.alpha*kx
    c_y=2*p.alpha*yphys
    dx_x=np.zeros_like(kx)
    dx_y=-np.ones_like(kx)
    dy_x=np.ones_like(kx)
    dy_y=np.zeros_like(kx)
    R=np.sqrt(dx*dx+dy*dy+dz*dz)
    assert np.min(R)>0, 'grid contains band touching; use offset/nonzero quadrature'
    Ax=dx*dx_x+dy*dy_y+dz*dz_x
    Ay=dx*dx_y+dy*dy_y+dz*dz_y
    energy=np.stack((c-R,c+R),axis=0)
    vx=np.stack((c_x-Ax/R,c_x+Ax/R),axis=0)
    vy=np.stack((c_y-Ay/R,c_y+Ay/R),axis=0)
    if not shifted:
        dbx=-np.ones_like(kx); dby=np.zeros_like(kx); dbz=np.zeros_like(kx)
        cb=np.zeros_like(kx); cxb=np.zeros_like(kx);cyb=np.zeros_like(kx)
        dxb_x=np.zeros_like(kx); dxb_y=np.zeros_like(kx)
        dzb_x=np.zeros_like(kx); dzb_y=np.zeros_like(kx)
    else:
        dbx=np.zeros_like(kx); dby=np.zeros_like(kx)
        dbz=dzb
        cb=-2*p.alpha*ky
        cxb=np.zeros_like(kx); cyb=-2*p.alpha*np.ones_like(kx)
        dxb_x=np.zeros_like(kx);dxb_y=np.zeros_like(kx)
        dzb_x=dzxb
        dzb_y=dzyb
    DR=dx*dbx+dy*dby+dz*dbz
    Dx=dbx*dx_x+dby*dy_x+dbz*dz_x+dx*dxb_x+dz*dzb_x
    Dy=dbx*dx_y+dby*dy_y+dbz*dz_y+dx*dxb_y+dz*dzb_y
    eb=np.stack((cb-DR/R,cb+DR/R),axis=0)
    vxb=np.stack((cxb-Dx/R+Ax*DR/R**3,cxb+Dx/R-Ax*DR/R**3),axis=0)
    vyb=np.stack((cyb-Dy/R+Ay*DR/R**3,cyb+Dy/R-Ay*DR/R**3),axis=0)
    return energy,vx,vy,eb,vxb,vyb


def occupancy_weight(energy,p):
    x=(energy-p.mu)/(2*p.temp)
    q=np.exp(-2*np.abs(x))
    w=q/(p.temp*(1+q)**2)
    wp=-w*np.tanh(x)/p.temp
    return w,wp


def full_response(p:Surface,b=0.,*,nrad=140,nang=256,shifted=False):
    kx,ky,measure=grid(nrad,nang,p.radius)
    e,vx,vy,*_=_elements(p,kx,ky,b=b,shifted=shifted)
    w,_=occupancy_weight(e,p)
    return float(np.sum(measure*vx*vy*w))


def analytic_terms(p:Surface,*,nrad=140,nang=256,shifted=False):
    kx,ky,measure=grid(nrad,nang,p.radius)
    e,vx,vy,eb,vxb,vyb=_elements(p,kx,ky,b=0.,shifted=shifted)
    w,wp=occupancy_weight(e,p)
    d={'velocity_x':float(np.sum(measure*vxb*vy*w)),
       'velocity_y':float(np.sum(measure*vx*vyb*w)),
       'population':float(np.sum(measure*vx*vy*wp*eb))}
    d['sum']=sum(d.values())
    return d


def finite_slope(p:Surface,h=.001,*,nrad=140,nang=256,shifted=False):
    if h<=0: raise ValueError('positive difference step required')
    a=full_response(p,h,nrad=nrad,nang=nang,shifted=shifted)
    b=full_response(p,-h,nrad=nrad,nang=nang,shifted=shifted)
    return (a-b)/(2*h)


def shift_vertex_residual(p:Surface,kx,ky):
    scalar=-2*p.alpha*ky
    z=(-3*p.lam*(ky*ky-kx*kx) if p.warp_axis=="mirror_y" else 6*p.lam*kx*ky)
    return scalar,z


def matrix_h(p:Surface,kx,ky,b=0.):
    kx=np.asarray(kx);ky=np.asarray(ky)
    c=p.beta*kx+p.alpha*(kx*kx+ky*ky)
    dz=(p.lam*(ky**3-3*ky*kx*kx) if p.warp_axis=="mirror_y" else p.lam*(kx**3-3*kx*ky*ky))
    dx=-(ky+b);dy=kx
    return np.stack([np.stack([c+dz,dx-1j*dy],axis=-1),
                     np.stack([dx+1j*dy,c-dz],axis=-1)],axis=-2)


def matrix_vy(p:Surface,kx,ky):
    c=2*p.alpha*ky
    dz=(3*p.lam*(ky*ky-kx*kx) if p.warp_axis=="mirror_y" else -6*p.lam*kx*ky)
    return np.stack([np.stack([c+dz,-np.ones_like(kx)],axis=-1),
                     np.stack([-np.ones_like(kx),c-dz]],axis=-1)],axis=-2)


def matrix_zeeman_b(p:Surface,kx,ky):
    kx=np.asarray(kx);ky=np.asarray(ky)
    return np.stack([np.stack([np.zeros_like(kx),-np.ones_like(kx)],axis=-1),
                     np.stack([-np.ones_like(kx),np.zeros_like(kx)],axis=-2)],axis=-2)
