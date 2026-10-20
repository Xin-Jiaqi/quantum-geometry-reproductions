"""Periodic lattice Pauli d-vector model extracted from the frozen v30 multiband model."""
from dataclasses import dataclass
import numpy as np
from ._spin_grid import I,SX,SY,SZ,mesh,fermi

@dataclass(frozen=True)
class LatticeModel:
    mass:float=.70
    tilt_sin_x:float=.23
    tilt_cos_y:float=.12
    xy_mix:float=.16
    y_cos_x:float=.08
    spin_orbit:float=1.

def fields(kx,ky,model=LatticeModel(),h=0.):
    x,y=kx,ky;p=model
    d=np.stack([p.spin_orbit*(np.sin(x)+p.xy_mix*np.sin(y)),
                p.spin_orbit*(np.sin(y)+p.y_cos_x*np.cos(x)),
                p.mass+np.cos(x)+np.cos(y)-h])
    dx=np.stack([p.spin_orbit*np.cos(x),-p.spin_orbit*p.y_cos_x*np.sin(x),-np.sin(x)])
    dy=np.stack([p.spin_orbit*p.xy_mix*np.cos(y),p.spin_orbit*np.cos(y),-np.sin(y)])
    dxy=np.zeros_like(d)
    dyy=np.stack([-p.spin_orbit*p.xy_mix*np.sin(y),-p.spin_orbit*np.sin(y),-np.cos(y)])
    e0=p.tilt_sin_x*np.sin(x)+p.tilt_cos_y*np.cos(y)
    e0x=p.tilt_sin_x*np.cos(x);e0y=-p.tilt_cos_y*np.sin(y)
    e0xy=np.zeros_like(x);e0yy=-p.tilt_cos_y*np.cos(y)
    return d,dx,dy,dxy,dyy,e0,e0x,e0y,e0xy,e0yy
