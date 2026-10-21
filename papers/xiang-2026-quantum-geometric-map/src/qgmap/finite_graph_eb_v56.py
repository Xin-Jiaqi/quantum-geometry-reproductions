"""v56 independent finite-graph E/B linear-response and mixed derivative.

This module does NOT calculate an intrinsic bulk orbital MNHE. It computes the
fully defined electric response of a finite, open Peierls graph with a scalar
site polarization, equilibrium f(H_B), and a single phenomenological coherence
rate gamma, in units hbar=1 and signed q=q_over_hbar.

rho_E = (gamma-i*omega+i[H_B,.])^-1 (i[P_b,rho_eq(H_B)]);
J_a = i[H_B,P_a];  sigma_ab=Tr(J_a rho_E).
Differentiate in B including (i) the current/contact vertex, (ii) the
Fermi-equilibrium Frechet derivative and (iii) the magnetic Liouvillian.
The complete sum, not its individual gauge-dependent terms, is physical within
this specified bath protocol. Finite E field refers to length-gauge -E*P.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .embedding_magnetic_v55 import peierls_field


def fermi(e, mu, temp):
    if not np.isfinite([mu, temp]).all() or temp <= 0:
        raise ValueError('finite chemical potential and strictly positive temperature required')
    x=np.clip((np.asarray(e)-mu)/temp,-500,500)
    return 1/(1+np.exp(x))


def equilibrium_density(h,mu,temp):
    e,u=np.linalg.eigh(h)
    return (u*fermi(e,mu,temp)[None,:])@u.conj().T


def frechet_fermi(h,m,mu,temp):
    """Exact Frechet derivative of f(H) with respect to H along dH=M.

    Handles degeneracies with finite-temperature derivative, where any basis in
    the degenerate subspace gives the same matrix result.
    """
    e,u=np.linalg.eigh(h)
    f=fermi(e,mu,temp)
    fp=-f*(1-f)/temp
    dif=e[:,None]-e[None,:]
    ff=f[:,None]-f[None,:]
    divided=np.divide(ff,dif,out=np.zeros_like(dif),where=abs(dif)>1e-9)
    midpoint=(e[:,None]+e[None,:])/2
    fm=fermi(midpoint,mu,temp)
    fd=-fm*(1-fm)/temp
    divided=np.where(abs(dif)>1e-9,divided,fd)
    mb=u.conj().T@m@u
    return u@(divided*mb)@u.conj().T


def polarization(positions,charge,axis):
    x=np.asarray(positions,float)
    if x.ndim!=2 or x.shape[1]!=2 or len(x)<2 or not np.isfinite(x).all() or not np.isfinite(charge) or charge==0:
        raise ValueError('positions finite (N,2) and nonzero charge required')
    if axis not in ('x','y'):raise ValueError('axis x or y required')
    return np.diag(charge*x[:,0 if axis=='x' else 1]).astype(complex)


def comm(a,b):return a@b-b@a


def current(h,p):
    """True dP/dt, including full B-dependent hopping; hbar=1."""
    return 1j*comm(h,p)


def resolvent(h,source,z):
    """Spectral solution (z+i ad_H)^-1 source, z=gamma-i omega, gamma>0."""
    if np.real(z)<=0 or not np.isfinite([z.real,z.imag]).all():
        raise ValueError('gamma must be positive')
    e,u=np.linalg.eigh(h)
    ub=u.conj().T
    s=ub@source@u
    zden=z+1j*(e[:,None]-e[None,:])
    return u@(s/zden)@ub


@dataclass(frozen=True)
class Spec:
    mu:float=0.1
    temperature:float=0.4
    gamma:float=0.08
    omega:float=0.35
    q_over_hbar:float=-0.87
    field:float=0.27
    probe_axis:str='x'
    drive_axis:str='y'
    def __post_init__(self):
        vals=(self.mu,self.temperature,self.gamma,self.omega,self.q_over_hbar,self.field)
        if not np.isfinite(vals).all() or self.temperature<=0 or self.gamma<=0 or self.q_over_hbar==0:
            raise ValueError('finite nonzero charge and positive T/gamma required')
        if self.probe_axis not in ('x','y') or self.drive_axis not in ('x','y'):
            raise ValueError('axes x/y required')


def state(hoppings,positions,spec=Spec(),gauge='symmetric'):
    """Full finite B electric linear response: no truncation in B."""
    q=spec.q_over_hbar  # hbar=1 chosen by the model
    p=peierls_field(hoppings,positions,spec.field,q,gauge=gauge)
    h,m=p.hamiltonian,p.magnetic_vertex
    pa=polarization(positions,q,spec.probe_axis)
    pb=polarization(positions,q,spec.drive_axis)
    j=current(h,pa)
    jB=current(m,pa)  # mandatory magnetic current contact
    rho0=equilibrium_density(h,spec.mu,spec.temperature)
    rhoB=frechet_fermi(h,m,spec.mu,spec.temperature)
    source=1j*comm(pb,rho0)
    z=spec.gamma-1j*spec.omega
    rhoE=resolvent(h,source,z)
    sig=np.trace(j@rhoE)
    return {'H':h,'M':m,'P_current':pa,'P_drive':pb,'J':j,'J_B':jB,
            'rho0':rho0,'rho_B':rhoB,'source_E':source,'rho_E':rhoE,
            'sigma':complex(sig),'z':z,'acceptance':'CONDITIONAL_FINITE_OPEN_GRAPH_SCALAR_GAMMA'}


def independent_sigma_spectral(h,pa,pb,mu,temp,gamma,omega):
    """Independent explicit double-band sum, no resolvent routine."""
    e,u=np.linalg.eigh(h)
    occupations=fermi(e,mu,temp)
    jb=u.conj().T@current(h,pa)@u
    pp=u.conj().T@pb@u
    total=0j
    for m in range(len(e)):
        for n in range(len(e)):
            source=1j*pp[m,n]*(occupations[n]-occupations[m])
            total+=jb[n,m]*source/(gamma-1j*omega+1j*(e[m]-e[n]))
    return complex(total)


def mixed_eb(hoppings,positions,spec=Spec(),gauge='symmetric'):
    """Exact B derivative of finite AC electric response incl. three terms.

    In a magnetic gauge G(B), current-contact, f(H)' and Liouville terms
    separately transform with derivative connections. Only their sum is gauge
    invariant; this is NOT a standalone intrinsic bulk magnetoelectric tensor.
    """
    s=state(hoppings,positions,spec,gauge)
    h,m,pa,pb=s['H'],s['M'],s['P_current'],s['P_drive']
    z=s['z']
    rhoE=s['rho_E']
    j=s['J']
    drho_source=resolvent(h,1j*comm(pb,s['rho_B']),z)
    drho_dynamics=resolvent(h,-1j*comm(m,rhoE),z)
    contributions={
        'current_contact':complex(np.trace(s['J_B']@rhoE)),
        'equilibrium_fermi_derivative':complex(np.trace(j@drho_source)),
        'magnetic_liouvillian':complex(np.trace(j@drho_dynamics))}
    combined=sum(contributions.values())
    return {'response_sigma':s['sigma'],'mixed_dsigma_dB':complex(combined),
            'contributions':contributions,
            'acceptance':'CONDITIONAL_FINITE_OPEN_GRAPH_SCALAR_GAMMA_MIXED_EB',
            'limitation':'Not bulk MNHE: finite OPEN graph, finite-frequency scalar gamma bath; components not individually gauge-invariant'}


def finite_field_derivative(hoppings,positions,spec=Spec(),gauge='symmetric',step=1e-4,points=5):
    """Independent full finite-B conductivity FD; 3/5-point stencils."""
    if not np.isfinite(step) or step<=0 or points not in (3,5):raise ValueError('step>0 and 3/5 stencil')
    from dataclasses import replace
    def calc(b):return state(hoppings,positions,replace(spec,field=b),gauge)['sigma']
    b=spec.field
    if points==3:return (calc(b+step)-calc(b-step))/(2*step)
    return (calc(b-2*step)-8*calc(b-step)+8*calc(b+step)-calc(b+2*step))/(12*step)


def equilibrium_current(hoppings,positions,spec=Spec(),gauge='symmetric',E=0.0):
    """Exact stationary equilibrium current for any scalar static E and B.

    J(B)=i[H(B),P]; after scalar E adds -E P_b, [P_a,P_b]=0,
    therefore Tr(f(H(B)-E P_b) J(B))=0 in every finite open graph.
    """
    from dataclasses import replace
    p=state(hoppings,positions,spec,gauge)
    field_h=p['H']-E*p['P_drive']
    rho=equilibrium_density(field_h,spec.mu,spec.temperature)
    return complex(np.trace(p['J']@rho))
