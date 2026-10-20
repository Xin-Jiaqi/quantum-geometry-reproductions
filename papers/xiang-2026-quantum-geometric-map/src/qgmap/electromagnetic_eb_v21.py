"""v21: finite-size exact spin-B derivative of electromagnetic *linear* Kubo response.

H(B)=H0+B Z, with site-local Zeeman Z=-sigma_x. The electric probe is
an explicitly finite-wavevector scalar potential phi(t)*rho_{-q}. For fixed
real omega+i delta, d_B chi_{A,rho_-q} is a mixed (phi,B) coefficient.

Two independent methods: (i) finite B derivatives of exact equilibrium
Lehmann susceptibility; (ii) linearized Liouville resolvent including
both [V, rho_eq'] and [Z,delta_rho_E]. The equilibrium derivative
rho_eq'=Df(H)[Z] keeps diagonal Fermi-surface f' terms. Freezing occupations
is a NAMED alternative and does not reproduce the equilibrated result.

Ward contact is independently differentiated for current source. Explicit
pure-gauge Peierls H(a)=e^{iaF} H e^{-iaF} includes current contact term;
zero physical static longitudinal response is an exact finite-size test.
Not Xiang's homogeneous DC transverse PHE or microscopic impurity scattering.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .real_frequency_continuity_v20 import ExactSystem,exact_hamiltonian,charge_operators,fermi,exact_retarded_chi,equilibrium_contact

@dataclass(frozen=True)
class EBSpec:
    system:ExactSystem=ExactSystem(L=4)
    direction:int=1
    omega:float=.31
    delta:float=.045
    magnetic_step:float=1e-3
    gauge_step:float=2e-4
    def __post_init__(self):
        if self.direction not in (0,1) or self.delta<=0 or self.magnetic_step<=0 or self.gauge_step<=0:
            raise ValueError('invalid response input')


def zeeman_operator(p:ExactSystem):
    """Uniform physical Pauli sigma_x; B measured in energy units, not Tesla."""
    sx=np.array([[0,1],[1,0]],complex)
    return -np.kron(np.eye(p.L*p.L),sx)


def build(spec:EBSpec):
    h,pos,_=exact_hamiltonian(spec.system)
    Z=zeeman_operator(spec.system)
    rho,rhobar,j,jbar,Q=charge_operators(h,pos,direction=spec.direction)
    # onsite spin Zeeman commutes with local charge density: J(B)=J(0)
    if np.max(abs(Z@rho-rho@Z))>1e-13: raise AssertionError('Zeeman must commute with rho_q')
    return h,Z,rho,rhobar,j,jbar,Q,pos


def equilibrium_response(h,A,V,p,omega,delta):
    e,u=np.linalg.eigh(h)
    f=fermi(e,p.mu,p.temperature)
    return complex(exact_retarded_chi(e,u,A,V,f,omega,delta))


def finite_zeeman_slope(spec:EBSpec,Aname='current',Vname='density',occupation_protocol='reequilibrated',step=None):
    """Central B derivative; frozen mode keeps B=0 eigenstate occupations ONLY
    as a diagnostic (not a gauge-covariant interacting preparation protocol)."""
    h,Z,rho,rhobar,j,jbar,Q,pos=build(spec)
    A={'density':rho,'current':j}[Aname]
    V={'density':rhobar,'current':jbar}[Vname]
    st=step or spec.magnetic_step
    if occupation_protocol=='reequilibrated':
        xp=equilibrium_response(h+st*Z,A,V,spec.system,spec.omega,spec.delta)
        xm=equilibrium_response(h-st*Z,A,V,spec.system,spec.omega,spec.delta)
    elif occupation_protocol=='frozen_label':
        ef,_=np.linalg.eigh(h)
        f0=fermi(ef,spec.system.mu,spec.system.temperature)
        def at(sign):
            es,us=np.linalg.eigh(h+sign*st*Z)
            return exact_retarded_chi(es,us,A,V,f0,spec.omega,spec.delta)
        xp,xm=at(1),at(-1)
    else:
        raise ValueError('unknown protocol')
    return (xp-xm)/(2*st)


def derivative_equilibrium_density(e,u,z,mu,T):
    """Fréchet derivative D f(H)[Z] including n=m term f'(energy)."""
    n=len(e); f=fermi(e,mu,T)
    en=e[:,None]-e[None,:]
    fn=f[:,None]-f[None,:]
    quot=np.empty((n,n),float)
    mask=np.abs(en)>1e-9
    np.divide(fn,en,out=quot,where=mask)
    # continuously define n=m and exact degeneracies by derivative
    mid=(e[:,None]+e[None,:])*.5
    fmid=fermi(mid,mu,T)
    quot[~mask]=(-fmid*(1-fmid)/T)[~mask]
    zb=u.conj().T@z@u
    return quot*zb


def liouville_slope(spec:EBSpec,Aname='current',Vname='density',include_equilibrium_derivative=True):
    h,Z,rho,rhobar,j,jbar,Q,pos=build(spec)
    A={'density':rho,'current':j}[Aname]
    V={'density':rhobar,'current':jbar}[Vname]
    e,u=np.linalg.eigh(h)
    f=fermi(e,spec.system.mu,spec.system.temperature)
    rho0=np.diag(f)
    zb=u.conj().T@Z@u
    vb=u.conj().T@V@u
    ab=u.conj().T@A@u
    den=spec.omega+1j*spec.delta-(e[:,None]-e[None,:])
    comm=lambda a,b:a@b-b@a
    x0=comm(vb,rho0)/den
    if include_equilibrium_derivative:
        drho=derivative_equilibrium_density(e,u,Z,spec.system.mu,spec.system.temperature)
    else:
        drho=np.zeros_like(zb)
    x1=(comm(vb,drho)+comm(zb,x0))/den
    return complex(np.trace(ab@x1))


def mixed_ward(spec:EBSpec,protocol='liouville'):
    h,Z,rho,rhobar,j,jbar,Q,pos=build(spec)
    if protocol=='liouville':
        calc=lambda a,b:liouville_slope(spec,a,b)
        # equilibrium magnetic derivative of equal-time contact
        e,u=np.linalg.eigh(h)
        d=u@derivative_equilibrium_density(e,u,Z,spec.system.mu,spec.system.temperature)@u.conj().T
        contact=np.trace(d@(rho@jbar-jbar@rho))
    elif protocol=='finite_difference':
        calc=lambda a,b:finite_zeeman_slope(spec,a,b)
        s=spec.magnetic_step
        def c(b):
            e,u=np.linalg.eigh(h+b*Z)
            f=fermi(e,spec.system.mu,spec.system.temperature)
            return equilibrium_contact(u,f,rho,jbar)
        contact=(c(s)-c(-s))/(2*s)
    else:raise ValueError('unknown protocol')
    w=spec.omega+1j*spec.delta
    density_res=w*calc('density','density')-Q*calc('current','density')
    current_res=w*calc('density','current')-Q*calc('current','current')-contact
    return {'density_source_ward_abs':float(abs(density_res)),
            'current_source_ward_abs':float(abs(current_res)),
            'contact_magnetic_derivative':complex(contact),
            'density_mixed':calc('density','density'),
            'current_mixed':calc('current','density')}


def gauge_contact(spec:EBSpec,b=0.,direction=None):
    """Static longitudinal pure-gauge Peierls A=grad F.

    Exact H(a)=U(a)H(B)U(a)^dagger, J(a)=U(a)J U(a)^dagger.
    Paramagnetic and explicit current/contact derivatives cancel.
    """
    h,Z,rho,rhobar,j,jbar,Q,pos=build(spec)
    f=np.cos((2*np.pi/spec.system.L)*pos[:,spec.direction])
    F=np.diag(np.repeat(f,2)).astype(complex)
    h=h+b*Z
    es,u=np.linalg.eigh(h)
    occ=fermi(es,spec.system.mu,spec.system.temperature)
    eq=(u*occ[None,:])@u.conj().T
    comm=lambda a,b:a@b-b@a
    deq=1j*comm(F,eq)
    dj=1j*comm(F,j)
    param=np.trace(deq@j)
    diamag=np.trace(eq@dj)
    a=spec.gauge_step
    theta=np.repeat(f,2)
    def current_at(amp):
        ph=np.exp(1j*amp*theta)
        gauged_h=ph[:,None]*h*ph.conj()[None,:]
        gauged_j=-(gauged_h@rho-rho@gauged_h)/Q
        vals,vecs=np.linalg.eigh(gauged_h)
        occb=fermi(vals,spec.system.mu,spec.system.temperature)
        return np.sum(occb*np.diag(vecs.conj().T@gauged_j@vecs))
    direct=(current_at(a)-current_at(-a))/(2*a)
    return {'paramagnetic':complex(param),'contact_diamagnetic':complex(diamag),
            'sum':complex(param+diamag),'full_finite_gauge_derivative':complex(direct)}


def gauge_mixed(spec:EBSpec):
    b=spec.magnetic_step
    plus=gauge_contact(spec,b)
    minus=gauge_contact(spec,-b)
    return {'paramagnetic_mixed':(plus['paramagnetic']-minus['paramagnetic'])/(2*b),
            'diamagnetic_mixed':(plus['contact_diamagnetic']-minus['contact_diamagnetic'])/(2*b),
            'total_mixed':(plus['sum']-minus['sum'])/(2*b),
            'finite_gauge_mixed':(plus['full_finite_gauge_derivative']-minus['full_finite_gauge_derivative'])/(2*b)}
