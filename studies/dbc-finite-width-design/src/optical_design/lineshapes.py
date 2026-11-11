"""Local quadratic-quartic response, line-shape dependent finite-width optimum.
Gaussian scale is its standard deviation; Lorentzian scale is its HWHM.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import pbdv
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad

@lru_cache(maxsize=32)
def gauss(order):return leggauss(order)

def profile_x(d,shape):
    d=np.asarray(d,float)
    if shape=='lorentzian':return np.real(1/np.sqrt(d+1j))
    if shape!='gaussian':raise ValueError('Unknown shape')
    out=np.zeros_like(d);mask=d>-38
    dv=d[mask];out[mask]=np.exp(-dv*dv/4)*pbdv(-.5,-dv)[0]/np.sqrt(2)
    return out

def phi(xi,zeta,shape='gaussian',order=280):
    if not np.isfinite([xi,zeta]).all():raise ValueError('Finite parameters needed')
    t,w=gauss(order)
    if shape=='gaussian':
        ymax=np.sqrt((zeta+np.sqrt(zeta*zeta+4*(16+abs(xi))))/2)
        y=(t+1)*ymax/2
        return float(ymax*np.sum(w*profile_x(zeta*y*y-y**4-xi,shape)))
    # Map [0,infinity) to [0,1); includes both signs of y.
    u=(t+1)/2;y=u/(1-u)
    return float(np.sum(w*profile_x(zeta*y*y-y**4-xi,shape)/(1-u)**2))

def phi_quad(xi,zeta,shape='gaussian'):
    """Independent adaptive 1D integration, including Lorentzian infinite tails."""
    def f(y):return float(profile_x(zeta*y*y-y**4-xi,shape))
    val,err=quad(f,0,np.inf,epsabs=1e-10,epsrel=1e-10,limit=350)
    return 2*val,2*err

def optimum(shape='gaussian',order=720):
    # Multi-start in both energy and unfolding coordinate, not imposed 2-sigma result.
    sols=[]
    for xi,zeta in [(-.1,2.8),(-1.,0.),(0.,5.),(-2.,1.)]:
        s=minimize(lambda z:-phi(z[0],z[1],shape,order),[xi,zeta],method='Nelder-Mead',
                   options={'xatol':2e-8,'fatol':2e-11,'maxiter':500})
        if not s.success:raise RuntimeError('Local shape optimization failed')
        sols.append(s)
    s=min(sols,key=lambda r:r.fun)
    crit=minimize_scalar(lambda xi:-phi(xi,0.,shape,order),bounds=(-6,1),method='bounded',options={'xatol':1e-10})
    widthfactor=2*np.sqrt(2*np.log(2)) if shape=='gaussian' else 2.
    return dict(shape=shape,xi=float(s.x[0]),zeta=float(s.x[1]),peak=float(-s.fun),
                critical_peak=float(-crit.fun),critical_xi=float(crit.x),
                peak_over_critical=float(s.fun/crit.fun),separation_over_scale=float(s.x[1]**2/4),
                separation_over_FWHM=float(s.x[1]**2/(4*widthfactor)),
                order=order,multistart_values=[float(-x.fun) for x in sols])
