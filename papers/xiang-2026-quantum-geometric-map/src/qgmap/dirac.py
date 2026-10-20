"""Tilted surface Dirac and Zeeman Berry curvature. Model energies eV; k Å^-1.

Pauli matrices are σ, not ℏσ/2. Eq(11) and literal SM Eq(63) are named
protocols, not silently merged. These are model kernels, NOT steady-state DC.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad

KB_EV_K=8.617333262145e-5
HBAR_EV_PS=6.582119569e-4
MU_B_EV_T=5.788381806e-5
I=np.eye(2, dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)

@dataclass(frozen=True)
class Parameters:
    beta:float=1/8
    vF_eVA:float=1.
    tau_ps:float=.5
    g_factor:float=10.
    Bx_T:float=1.
    def __post_init__(self):
        if not 0<self.beta<1 or self.vF_eVA<=0 or self.tau_ps<=0:
            raise ValueError("Require 0<beta<1 and positive velocity and relaxation time")

def _sign(helicity):
    if helicity=="main-eq10": return 1
    if helicity=="sm-eq62": return -1
    raise ValueError("helicity must be 'main-eq10' or 'sm-eq62'")

def hamiltonian(kx,ky,p=Parameters(),helicity="main-eq10"):
    s=_sign(helicity)
    return p.beta*p.vF_eVA*kx*I+s*p.vF_eVA*(kx*SY-ky*SX)

def energies(kx,ky,p=Parameters()):
    kx,ky=np.broadcast_arrays(np.asarray(kx,float),np.asarray(ky,float))
    r=np.hypot(kx,ky)
    center=p.beta*p.vF_eVA*kx
    return center-p.vF_eVA*r,center+p.vF_eVA*r

def local_zbc_analytic(kx,ky,helicity="main-eq10"):
    s=_sign(helicity)
    kx,ky=np.broadcast_arrays(np.asarray(kx,float),np.asarray(ky,float))
    r=np.hypot(kx,ky)
    with np.errstate(divide="ignore",invalid="ignore"):
        zx=-s*kx*ky/r**3
        zy=s*kx*kx/r**3
    return np.where(r>0,zx,np.nan),np.where(r>0,zy,np.nan)

def local_zbc_eigensolver(kx,ky,p=Parameters(),helicity="main-eq10",phases=(0.,0.)):
    """Independent matrix-element calculation, gauge invariant for U(1) phases."""
    if np.hypot(kx,ky)==0: raise ValueError("Nondegenerate bands needed; k=0 is a node")
    s=_sign(helicity)
    e,u=np.linalg.eigh(hamiltonian(kx,ky,p,helicity))
    u=u@np.diag(np.exp(1j*np.asarray(phases,float)))
    sx=u.conj().T@SX@u
    vx=u.conj().T@(p.beta*p.vF_eVA*I+s*p.vF_eVA*SY)@u
    vy=u.conj().T@(-s*p.vF_eVA*SX)@u
    rx=1j*vx[0,1]/(e[1]-e[0])
    ry=1j*vy[0,1]/(e[1]-e[0])
    return float(-2*np.imag(rx*sx[1,0])),float(-2*np.imag(ry*sx[1,0]))

def eq11_kernel(beta):
    b=np.asarray(beta,float)
    if np.any((b<=0)|(b>=1)): raise ValueError("0 < beta < 1")
    return -b/(4*np.pi*(1+np.sqrt(1-b*b))**2)

def eq63_literal_integral(beta):
    if not 0<beta<1: raise ValueError("0 < beta < 1")
    val,err=quad(lambda theta: np.cos(theta)*np.sin(theta)**2/(1+beta*np.cos(theta)),
                 0,2*np.pi,epsabs=2e-13,epsrel=2e-13,limit=250)
    return float(-val/(4*np.pi**2)),float(err/(4*np.pi**2))

def continuum_temperature_kernel(mu_eV,temperature_K,beta,protocol="paper-eq11"):
    if protocol=="paper-eq11": k=eq11_kernel(beta)
    elif protocol=="sm-eq63-literal": k=eq63_literal_integral(beta)[0]
    else: raise ValueError("Unrecognized printed formula")
    mu=np.asarray(mu_eV,float)
    if temperature_K<0: raise ValueError("T must be nonnegative")
    if temperature_K==0:
        if np.any(mu==0):raise ValueError("Undefined T=0 mu=0 Dirac step")
        return k*np.sign(mu)
    return k*np.tanh(mu/(2*KB_EV_K*temperature_K))

def zeeman_time_prefactor(p=Parameters()):
    """g μ_B B τ/ℏ, using angular frequency, not cyclic frequency."""
    return p.g_factor*MU_B_EV_T*p.Bx_T*p.tau_ps/HBAR_EV_PS

def equilibrium_drude_shift_invariant(p=Parameters(),mu_eV=.012,T_K=20.,n=121,cutoff_Ainv=.15):
    """Odd-q_y full-plane, shifted-momentum reference control for uniform B_x."""
    if T_K<=0 or n<3: raise ValueError("Requires finite T and n >= 3")
    grid=np.linspace(-cutoff_Ainv,cutoff_Ainv,n)
    x,y=np.meshgrid(grid,grid,indexing="ij")
    k=np.hypot(x,y);ksafe=np.where(k==0,1,k)
    out=np.zeros_like(x)
    for s in (-1,1):
        e=p.beta*p.vF_eVA*x+s*p.vF_eVA*k
        vx=p.vF_eVA*(p.beta+s*x/ksafe)
        vy=p.vF_eVA*s*y/ksafe
        q=np.clip((e-mu_eV)/(2*KB_EV_K*T_K),-250,250)
        out+=vx*vy/(4*KB_EV_K*T_K*np.cosh(q)**2)
    return float(out.sum()*(grid[1]-grid[0])**2/(2*np.pi)**2)
