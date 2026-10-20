"""Direct SI-units EB Liouville and Eq(20) current for the paper's tilted Dirac model.

H and all resolvent denominators are in J, k in inverse METERS; derivatives dH/dk
in J m, electric potential +e E y in J, Zeeman H_B=-g mu_B B sigma_x
in J. Charge-current operator is +e/hbar dH/dk, for authors' e>0
convention; changing the electron charge convention flips current signs but NOT
magnitudes or the 2*pi question. Wave-vector measure is d^2k/(2pi)^2.

The band-resolved dynamic response below uses nonzero energy broadenings
eta_{E,B}=hbar/tau. No collision model is inferred from this prescription.
Not a claim of the measurable steady-state dc response.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import constants as c
from .si_fig2_v25 import geometry

I=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
JEV=c.e # J per eV, exact by SI definition

@dataclass(frozen=True)
class SIConfig:
    beta:float=1/8
    vF_eVnm:float=.30
    mu_eV:float=.060
    thermal_eV:float=.006
    kmax_eV:float=.48
    tau_s:float=5e-13
    E_Vm:float=1e4
    B_T:float=1.
    g:float=10.
    nr:int=128
    nphi:int=160
    eta_scale:float=1.
    omega_e_scale:float=0.
    omega_b_scale:float=0.

    def __post_init__(self):
        if not abs(self.beta)<1 or min(self.vF_eVnm,self.thermal_eV,self.kmax_eV,self.tau_s,self.nr,self.nphi,self.eta_scale)<=0:
            raise ValueError('invalid model parameters')

def kgrid(p:SIConfig):
    vf=p.vF_eVnm*JEV*1e-9 # J*m, includes hbar implicitly
    cutoff=p.kmax_eV*JEV/vf
    x,w=leggauss(p.nr)
    r=cutoff*(x+1)/2
    rw=cutoff*w/2
    phi=(np.arange(p.nphi)+.317)*2*np.pi/p.nphi
    kx=(r[:,None]*np.cos(phi)[None,:]).ravel()
    ky=(r[:,None]*np.sin(phi)[None,:]).ravel()
    weights=np.broadcast_to((rw*r)[:,None]/(2*np.pi*p.nphi),(p.nr,p.nphi)).ravel()
    return np.stack([kx,ky],axis=-1),weights,vf

def _comm(a,b):return a@b-b@a

def _bands(p:SIConfig,k,vf):
    kx,ky=k.T
    H=vf*(p.beta*kx[:,None,None]*I+kx[:,None,None]*SY-ky[:,None,None]*SX)
    vx=vf*(p.beta*I+SY)
    vy=-vf*SX
    e,u=np.linalg.eigh(H)
    ud=u.conj().transpose(0,2,1)
    f=1/(1+np.exp(np.clip((e-p.mu_eV*JEV)/(p.thermal_eV*JEV),-500,500)))
    fp=-f*(1-f)/(p.thermal_eV*JEV)
    delta=e[:,:,None]-e[:,None,:]
    ratio=np.zeros_like(delta)
    off=np.abs(delta)>1e-26
    num=f[:,:,None]-f[:,None,:]
    ratio[off]=num[off]/delta[off]
    ratio[:,0,0]=fp[:,0];ratio[:,1,1]=fp[:,1]
    vy_band=ud@vy@u
    drho=u@(ratio*vy_band)@ud # m, derivative w.r.t physical k (1/m)
    rho=u@(np.einsum('ni,ij->nij',f,I))@ud
    return e,u,ud,rho,drho,vx,vy,vy_band,f,fp

def _sum_current(rho,vx,weights):
    # sheet current density A/m
    return c.e/c.hbar*complex(np.sum(weights*np.trace(vx@rho,axis1=1,axis2=2)))

def eq20_physical(p:SIConfig):
    """Literal SM Eq20 band sum with all constants restored independently of F(beta)."""
    k,wt,vf=kgrid(p)
    e,u,ud,rho,drho,vx,vy,vy_band,f,fp=_bands(p,k,vf)
    z=-p.g*c.physical_constants['Bohr magneton'][0]*p.B_T*SX
    zb=ud@z@u
    vxb=ud@vx@u
    df=fp*np.real(np.diagonal(vy_band,axis1=1,axis2=2))
    dfnm=df[:,None,:]-df[:,:,None]
    ediff=e[:,:,None]-e[:,None,:]
    # eps_{mn}: careful: in this routine epsilon_{nm}=eps_n-eps_m,
    # follows the v24 direct rho matrix convention.
    den=np.where(np.eye(2,dtype=bool)[None,:,:],1.,ediff)
    rho_band=-(c.e*p.E_Vm*p.tau_s/c.hbar)*zb*dfnm/den
    rho_band[:,0,0]=0;rho_band[:,1,1]=0
    return c.e/c.hbar*complex(np.sum(wt*np.einsum('nij,nji->n',rho_band,vxb)))

def two_color_physical(p:SIConfig):
    """Both field orderings from resolvents i*hbar*d_rho=[H,rho]."""
    k,wt,vf=kgrid(p)
    e,u,ud,rho,drho,vx,vy,_,_,_=_bands(p,k,vf)
    eb=e[:,:,None]-e[:,None,:]
    z=-p.g*c.physical_constants['Bohr magneton'][0]*p.B_T*SX
    eta=c.hbar/p.tau_s*p.eta_scale
    omE=p.omega_e_scale*eta;omB=p.omega_b_scale*eta
    def R(X,omega,etaE):
        return u@((ud@X@u)/(omega+1j*etaE-eb))@ud
    RE=R(1j*c.e*p.E_Vm*drho,omE,eta)
    RB=R(_comm(z,rho),omB,eta)
    dRB=R(_comm(z,drho)+_comm(vy,RB),omB,eta)
    BE=R(1j*c.e*p.E_Vm*dRB,omE+omB,2*eta)
    EB=R(_comm(z,RE),omE+omB,2*eta)
    return {'b_then_e_A_per_m':_sum_current(BE,vx,wt),
            'e_then_b_A_per_m':_sum_current(EB,vx,wt),
            'sum_A_per_m':_sum_current(BE+EB,vx,wt),
            'eta_energy_J':eta,'F_exact':-p.beta/(4*np.pi*(1+np.sqrt(1-p.beta**2))**2),
            'natural_dimensionless_tau_Zeeman_phase':p.tau_s*p.g*c.physical_constants['Bohr magneton'][0]*p.B_T/c.hbar}

def eq11_SI_magnitude(p:SIConfig):
    T=float(np.tanh(p.mu_eV/(2*p.thermal_eV)))
    return abs(geometry(p.beta)*T)*p.tau_s*c.e**2/c.hbar**2*p.g*c.physical_constants['Bohr magneton'][0]*abs(p.B_T)

def frequency_hypotheses(p:SIConfig):
    # Conventions never interchangeable unless a physical derivation proves it.
    ez=p.g*c.physical_constants['Bohr magneton'][0]*p.B_T
    return {'E_Z_J':ez, 'omega_Z_rad_s':ez/c.hbar,'f_Z_Hz':ez/c.h,
            'tau_omega_Z_rad':p.tau_s*ez/c.hbar,
            'tau_f_Z_cycles':p.tau_s*ez/c.h,
            'rad_per_cycle':(ez/c.hbar)/(ez/c.h)}
