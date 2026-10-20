"""Continuum tilted TI surface, homogeneous length-gauge two-color EB response.

Reproduce source SM Eqs. (4)-(10),(20) at fixed continuum UV cutoff using
independent spectral-projector and direct eigensystem algorithms.

H_main = beta kx I + kx sigma_y - ky sigma_x (main Eq. 10).
H_sm   = beta kx I + ky sigma_x - kx sigma_y (SM Eq. 62).
A magnetic 'field' bbar multiplies sigma_x; physical B has bbar=-g mu_B B.
e=hbar=v_F=1; weights dkx dky/(2pi)^2; the complex E_y and bbar
amplitudes are 1, without the original real-field 1/2 frequency factors.
No relaxation collision operator or physically converged DC transport claim.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
I=np.eye(2, dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)

@dataclass(frozen=True)
class DiracSpec:
    beta:float=.125
    mu:float=.20
    temp:float=.035
    cutoff:float=2.5
    nr:int=96
    nphi:int=128
    eta_e:float=.012
    eta_b:float=.012
    omega_e:float=0.
    omega_b:float=0.
    helicity:str='main'
    def __post_init__(self):
        if not abs(self.beta)<1 or not self.temp>0 or self.cutoff<=0 or self.nr<6 or self.nphi<16:
            raise ValueError('Require |beta|<1, T>0, cutoff>0, nr>=6, nphi>=16')
        if self.eta_e<=0 or self.eta_b<=0 or self.helicity not in ('main','sm'):
            raise ValueError('positive eta and helicity main/sm required')

@lru_cache(maxsize=32)
def quadrature(nr,nphi,cutoff,angular_shift=.37):
    x,w=leggauss(nr)
    r=(x+1)*cutoff/2
    rw=w*cutoff/2
    phi=(np.arange(nphi)+angular_shift)*2*np.pi/nphi
    kx=(r[:,None]*np.cos(phi)[None,:]).reshape(-1)
    ky=(r[:,None]*np.sin(phi)[None,:]).reshape(-1)
    weights=np.broadcast_to((rw*r)[:,None]*(1/(2*np.pi*nphi)),(nr,nphi)).ravel()
    return np.stack([kx,ky],axis=-1),weights

def cb(x,y):return x@y-y@x

def fermi(e,mu,t):return 1/(1+np.exp(np.clip((e-mu)/t,-500,500)))

def traces(weight,vx,rho):return complex(np.sum(weight*np.trace(vx@rho,axis1=-2,axis2=-1)))

def continuum_matrices(spec,k):
    x,y=k.T
    sign=1 if spec.helicity=='main' else -1
    h=spec.beta*x[:,None,None]*I+sign*(x[:,None,None]*SY-y[:,None,None]*SX)
    vx=spec.beta*I+sign*SY
    vy=-sign*SX
    return h,vx,vy,SX

def analytic_projectors(spec,k):
    x,y=k.T;r=np.hypot(x,y)
    sign=1 if spec.helicity=='main' else -1
    nmat=sign*(x[:,None,None]*SY-y[:,None,None]*SX)/r[:,None,None]
    dn=-sign*SX/r[:,None,None]-nmat*(y/(r*r))[:,None,None]
    pp=(I+nmat)*.5;pm=(I-nmat)*.5
    dp=dn*.5
    e=np.stack((spec.beta*x+r,spec.beta*x-r),axis=1)
    projectors=np.stack((pp,pm),axis=1)
    f=fermi(e,spec.mu,spec.temp)
    fp=-f*(1-f)/spec.temp
    velocity_y=np.stack((y/r,-y/r),axis=1)
    rho0=np.einsum('na,naij->nij',f,projectors)
    drho0=np.einsum('na,naij->nij',fp*velocity_y,projectors)+(f[:,0]-f[:,1])[:,None,None]*dp
    drho_pop=np.einsum('na,naij->nij',fp*velocity_y,projectors)
    return e,projectors,rho0,drho0,drho_pop

def projector_resolvent(source,e,projectors,omega,eta):
    # band eigenvalues order (+,-); finite positive r by Gauss quadrature
    delta=e[:,:,None]-e[:,None,:]
    out=np.zeros_like(source)
    for a in range(2):
        for b in range(2):
            out+=projectors[:,a]@source@projectors[:,b] / (omega+1j*eta-delta[:,a,b])[:,None,None]
    return out

def projector_channels(spec:DiracSpec,detail=False):
    k,weight=quadrature(spec.nr,spec.nphi,spec.cutoff)
    h,vx,vy,z=continuum_matrices(spec,k)
    e,p,rho0,drho0,drho_pop=analytic_projectors(spec,k)
    RE=projector_resolvent(1j*drho0,e,p,spec.omega_e,spec.eta_e)
    RB=projector_resolvent(cb(z,rho0),e,p,spec.omega_b,spec.eta_b)
    # exact analytic derivative of magnetic resolvent; Z independent of k
    dRB=projector_resolvent(cb(z,drho0)+cb(vy,RB),e,p,spec.omega_b,spec.eta_b)
    BE=projector_resolvent(1j*dRB,e,p,spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    EB=projector_resolvent(cb(z,RE),e,p,spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    RdiagE=projector_resolvent(1j*drho_pop,e,p,spec.omega_e,spec.eta_e)
    PoleChannel=projector_resolvent(cb(z,RdiagE),e,p,spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    result={"b_then_e":traces(weight,vx,BE),"e_then_b":traces(weight,vx,EB),
            "total":traces(weight,vx,BE+EB),"e_population_then_b":traces(weight,vx,PoleChannel)}
    if detail:result.update(k=k,weight=weight,BE=BE,EB=EB,drho0=drho0,e=e,p=p)
    return result

def eigen_channels(spec:DiracSpec):
    """Independent eigensystem, Fréchet derivative and band Green resolvents."""
    k,weight=quadrature(spec.nr,spec.nphi,spec.cutoff)
    h,vx,vy,z=continuum_matrices(spec,k)
    e,u=np.linalg.eigh(h)
    ud=np.swapaxes(u.conj(),-1,-2)
    eb=e[:,:,None]-e[:,None,:]
    f=fermi(e,spec.mu,spec.temp);fp=-f*(1-f)/spec.temp
    ratio=np.empty_like(eb)
    mask=np.abs(eb)>1e-12
    numerator=f[:,:,None]-f[:,None,:]
    ratio[mask]=numerator[mask]/eb[mask]
    ratio[:,0,0]=fp[:,0];ratio[:,1,1]=fp[:,1]
    vyb=ud@vy@u
    drho0=u@(ratio*vyb)@ud
    rho0=u@(f[:,:,None]*np.eye(2)[None,:,:])@ud
    def solve(s,om,eta):return u@((ud@s@u)/(om+1j*eta-eb))@ud
    re=solve(1j*drho0,spec.omega_e,spec.eta_e)
    rb=solve(cb(z,rho0),spec.omega_b,spec.eta_b)
    drb=solve(cb(z,drho0)+cb(vy,rb),spec.omega_b,spec.eta_b)
    be=solve(1j*drb,spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    er=solve(cb(z,re),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    return dict(b_then_e=traces(weight,vx,be),e_then_b=traces(weight,vx,er),total=traces(weight,vx,be+er))

def eq20_direct(spec:DiracSpec):
    """Current/tau, literal SM Eq.(20), field bbar (NOT physical B)."""
    k,weight=quadrature(spec.nr,spec.nphi,spec.cutoff)
    h,vx,vy,z=continuum_matrices(spec,k)
    e,u=np.linalg.eigh(h);ud=np.swapaxes(u.conj(),-1,-2)
    f=fermi(e,spec.mu,spec.temp)
    df=-f*(1-f)/spec.temp*np.real(np.diagonal(ud@vy@u,axis1=-2,axis2=-1))
    zband=ud@z@u;vband=ud@vx@u
    epsmn=e[:,:,None]-e[:,None,:]
    deriv_f_nm=df[:,None,:]-df[:,:,None]
    rho=-zband*deriv_f_nm/np.where(np.eye(2,dtype=bool)[None,:,:],1,epsmn)
    rho[:,0,0]=0;rho[:,1,1]=0
    return complex(np.sum(weight*np.einsum('nij,nji->n',rho,vband)))

def angular_exact(spec:DiracSpec):
    """Infinite-UV analytic radial integral of literal SM Eq20 for H_main.
    Due to helicity conventions, SM Eq 62 has the opposite sign of H_main
    in the bbar-normalized Zeeman-Berry coefficient.
    """
    beta=spec.beta
    # Eq20 direct: opposite of main Eq11 coefficient in the main-convention calculation
    if abs(beta)<1e-6:return 0.
    # Numerically stable expression: original numerator=(beta^2-2*(1-sqrt(1-beta^2)))
    # = - beta^4/(1+sqrt(1-beta^2))^2.
    closed=beta/(4*np.pi*(1+np.sqrt(1-beta**2))**2)
    return closed*np.tanh(spec.mu/(2*spec.temp))*(1 if spec.helicity=='main' else -1)
