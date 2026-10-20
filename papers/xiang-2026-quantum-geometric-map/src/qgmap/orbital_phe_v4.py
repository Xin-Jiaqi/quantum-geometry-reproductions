"""Orbital PHE, SM Eq. (53)-(56), two-band periodic 2D B_z subcase.

Scientific audit: compare Eq.56 first and second displayed lines independently with
finite-B semiclassical Boltzmann conductivity.  The three-band processes of SM
Eq.(49) are not evaluated. All numbers e=hbar=tau=1, k in radians, T in energy.
Model has an isolated two-band spectrum and a band-independent orbital moment
m_z=(epsilon_n-epsilon_other)*Omega_n/2.
"""
from __future__ import annotations
import numpy as np
from ._orbital_grid import LatticeModel, fields, mesh, fermi


def _dot(a,b):
    return np.einsum('i...,i...->...',a,b)


def _cross(a,b):
    return np.cross(a,b,axis=0)


def second_fermi(e,mu,T):
    f,fp=fermi(e,mu,T)
    return f,fp,f*(1-f)*(1-2*f)/T**2


def periodic_data(N=141,model=LatticeModel()):
    x,y=mesh(N)
    d,dx,dy,dxy,dyy,e0,e0x,e0y,e0xy,e0yy=fields(x,y,model)
    so=model.spin_orbit
    dxx=np.stack([-so*np.sin(x),-so*model.y_cos_x*np.cos(x),-np.cos(x)])
    K=_dot(d,_cross(dx,dy)); D=np.linalg.norm(d,axis=0)
    if float(np.min(D))<1e-8:raise ValueError('Degeneracy in band geometry: reject Abelian formulas')
    sign=np.stack([-np.ones_like(D),np.ones_like(D)],axis=-1)
    E=e0[...,None]+sign*D[...,None]
    gradx=_dot(d,dx);grady=_dot(d,dy)
    Vx=e0x[...,None]+sign*(gradx/D)[...,None]
    Vy=e0y[...,None]+sign*(grady/D)[...,None]
    Hxy=e0xy[...,None]+sign*((_dot(dx,dy)+_dot(d,dxy))/D-gradx*grady/D**3)[...,None]
    Hxx=(-model.tilt_sin_x*np.sin(x))[...,None]+sign*((_dot(dx,dx)+_dot(d,dxx))/D-gradx**2/D**3)[...,None]
    Hyy=e0yy[...,None]+sign*((_dot(dy,dy)+_dot(d,dyy))/D-grady**2/D**3)[...,None]
    # Berry curvature convention of main paper; H_B^orbital = -m_z B_z.
    om0=K/(2*D**3)
    Omega=np.stack([om0,-om0],axis=-1)
    moment=-K/(2*D**2)
    dKx=_dot(d,_cross(dxx,dy)+_cross(dx,dxy))
    dKy=_dot(d,_cross(dxy,dy)+_cross(dx,dyy))
    moment_x=-dKx/(2*D**2)+K*gradx/D**4
    moment_y=-dKy/(2*D**2)+K*grady/D**4
    return dict(E=E,v=[Vx,Vy],h=[[Hxx,Hxy],[Hxy,Hyy]],Omega=Omega,
                m=moment[...,None],dm=[moment_x[...,None],moment_y[...,None]],
                minimum_direct_gap=float(2*D.min()))


def orbital_phe_eq56(N=141,mu=1.15,T=.10,model=LatticeModel(),a=0,b=1,B_step=1e-3):
    """Keep SM Eq.56 first/second lines alongside independent IBP & finite-field results."""
    if a not in (0,1) or b not in (0,1):raise ValueError('2D in-plane directions only')
    if B_step<=0 or T<=0:raise ValueError('B_step and T must be positive')
    q=periodic_data(N,model)
    E=q['E']; vx=q['v'][a];vy=q['v'][b];h=q['h'][a][b]
    m=q['m'];dma=q['dm'][a];dmb=q['dm'][b];Omega=q['Omega']
    _,fp,fpp=second_fermi(E,mu,T)
    geom=-vx*vy*Omega*fp
    first= -dma*vy*fp -vx*(dmb*fp+m*vy*fpp) +geom
    second=2*m*h*fp -vx*vy*m*fpp+geom
    ibp=2*m*h*fp +vx*vy*m*fpp+geom
    def integrate(z):return float(np.mean(np.sum(z,axis=-1)))
    def sigma(B):
        eps=E-B*m
        _,deriv=fermi(eps,mu,T)
        return integrate((vx-B*dma)*(vy-B*dmb)*deriv/(1+B*Omega))
    fd=(sigma(B_step)-sigma(-B_step))/(2*B_step)
    return {'protocol':'orbital_PHE_SM_eq56_literal_vs_independent_finite_B_2D_Bz',
       'mesh':N,'mu':mu,'T':T,'B_step':B_step,'a':'xy'[a],'b':'xy'[b],
       'normalization':'sigma_ab,z^(1)/tau, e=hbar=1, BZ mean over [-pi,pi)^2',
       'minimum_direct_gap':q['minimum_direct_gap'],
       'eq56_first_line':integrate(first),'eq56_second_line_as_printed':integrate(second),
       'eq56_integrated_by_parts_corrected_sign':integrate(ibp),
       'independent_finite_B_boltzmann_derivative':fd,
       'fermi_second_derivative_integral':integrate(vx*vy*m*fpp),
       'Berry_curvature_part':integrate(geom),
       'printed_line_discrepancy':integrate(second-first),
       'ibp_minus_first':integrate(ibp-first),
       'finite_difference_minus_first':fd-integrate(first),
       'finite_difference_antisymmetry_step':B_step}
