"""Classical parabolic-band orbital OHE, NOT a multiband quantum-metric solver."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from .dirac import KB_EV_K

@dataclass(frozen=True)
class OHE:
    density:float
    sigma_xy_over_tau2:float
    sigma_yx_over_tau2:float
    relative_error:float

def ohe_parabolic_2d(mass=1.7,mu=.22,temperature_K=300.,cutoff_eV=None):
    if mass<=0 or temperature_K<=0: raise ValueError("Positive mass and temperature required")
    kbt=KB_EV_K*temperature_K
    cutoff=float(mu+45*kbt if cutoff_eV is None else cutoff_eV)
    if cutoff<=max(mu,0):raise ValueError("Cutoff too low")
    def f(e):
        return 1/(1+np.exp(np.clip((e-mu)/kbt,-500,500)))
    def weight(e):
        z=np.clip((e-mu)/(2*kbt),-300,300)
        return 1/(4*kbt*np.cosh(z)**2)
    density=mass/(2*np.pi)*quad(f,0,cutoff,epsabs=1e-13,epsrel=1e-13)[0]
    hall=quad(lambda e:e*weight(e),0,cutoff,epsabs=1e-13,epsrel=1e-13)[0]/(2*np.pi*mass)
    expectation=density/mass**2
    return OHE(float(density),float(hall),float(-hall),float(abs(hall-expectation)/abs(expectation)))
