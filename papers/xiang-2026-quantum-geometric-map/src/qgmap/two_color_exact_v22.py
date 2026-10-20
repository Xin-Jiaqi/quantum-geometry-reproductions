"""v22: two independently switched E and B drives in an exact finite system.

Two-color Liouville in the finite-q scalar-electric / uniform onsite Zeeman model:
  i dt rho = [H0 + E exp(s_E t) V + B exp(s_B t) Z, rho],
  s_a = eta_a - i omega_a, V=rho_-q, Z=-sigma_x.
Each complex drive amplitude is ONE (no 1/2 real-field Fourier convention).

At second order, rho_EB = R_sum([V,rho_B] + [Z,rho_E]);
R_a(X) = (omega_a+i eta_a-[H0,.])^(-1) X.
The two source terms are returned SEPARATELY.  The current operator is
J_q=-[H0,rho_q]/Q(q), independent of onsite B. The scalar V=rho_-q
is not the homogeneous length-gauge electric perturbation in Xiang SM (8).

Prepared-static-B susceptibility uses rho'_eq=D f(H0)[Z]; isolated,
adiabatically switched B has no diagonal occupation response.  Thus at
omega_B=0, eta_B -> 0, the difference is precisely the equilibrium f' Z_nn
population term, provided denominators remain regular at finite eta_E.

Not a microscopic scattering theory and NOT paper-defined DC PHE.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.linalg import solve_sylvester
from .real_frequency_continuity_v20 import ExactSystem, fermi
from .electromagnetic_eb_v21 import EBSpec, build, derivative_equilibrium_density


@dataclass(frozen=True)
class TwoColorSpec:
    system: ExactSystem = field(default_factory=lambda: ExactSystem(L=3))
    direction: int = 1
    omega_e: float = .23
    omega_b: float = .17
    eta_e: float = .14
    eta_b: float = .12

    def __post_init__(self):
        xs = (self.omega_e,self.omega_b,self.eta_e,self.eta_b)
        if self.direction not in (0,1) or not all(np.isfinite(x) for x in xs) or self.eta_e<=0 or self.eta_b<=0:
            raise ValueError('finite frequencies and strictly positive adiabatic regulators required')


def _data(spec: TwoColorSpec):
    # No performance shortcut: use exact finite disordered Hamiltonian.
    h,z,r,rm,j,jm,Q,pos=build(EBSpec(system=spec.system,direction=spec.direction))
    return h,z,rm,j


def _basis(spec):
    h,z,v,a=_data(spec)
    e,u=np.linalg.eigh(h)
    ud=u.conj().T
    f=fermi(e,spec.system.mu,spec.system.temperature)
    return e,ud@z@u,ud@v@u,ud@a@u,f,u


def _c(a,b): return a@b-b@a


def _r(e,x,omega,eta):
    return x/(omega+1j*eta-(e[:,None]-e[None,:]))


def frequency_channels(spec:TwoColorSpec):
    """Explicit E-on-B and B-on-E kernels in H0 eigenbasis."""
    e,z,v,a,f,u=_basis(spec)
    eq=np.diag(f.astype(complex))
    re=_r(e,_c(v,eq),spec.omega_e,spec.eta_e)
    rb=_r(e,_c(z,eq),spec.omega_b,spec.eta_b)
    # Here B->E denotes an E interaction acting on B-first density.
    b_then_e=_r(e,_c(v,rb),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    e_then_b=_r(e,_c(z,re),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    val=lambda x:complex(np.trace(a@x))
    return dict(b_then_e=val(b_then_e),e_then_b=val(e_then_b),total=val(b_then_e+e_then_b),
                rho_b=rb,rho_e=re,rho_eb=b_then_e+e_then_b,eig=e,zeeman=z,electric=v,observable=a,
                occupations=f)


def independent_sylvester_channels(spec:TwoColorSpec):
    """Independent orbital-basis linear-equation solution, not diagonal resolvent."""
    h,z,v,a=_data(spec)
    e,u=np.linalg.eigh(h)
    f=fermi(e,spec.system.mu,spec.system.temperature)
    rho0=(u*f[None,:])@u.conj().T
    eye=np.eye(len(h))
    def resolvent(src,w,eta):
        return solve_sylvester((w+1j*eta)*eye-h,h,src)
    re=resolvent(_c(v,rho0),spec.omega_e,spec.eta_e)
    rb=resolvent(_c(z,rho0),spec.omega_b,spec.eta_b)
    be=resolvent(_c(v,rb),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    eb=resolvent(_c(z,re),spec.omega_e+spec.omega_b,spec.eta_e+spec.eta_b)
    val=lambda x:complex(np.trace(a@x))
    return {'b_then_e':val(be),'e_then_b':val(eb),'total':val(be+eb)}


def static_preparation(spec:TwoColorSpec):
    """Exact equilibrium-B derivative of retarded E response, cf. v21."""
    if spec.omega_b!=0: raise ValueError('static magnetic preparation requires omega_b=0')
    e,z,v,a,f,u=_basis(spec)
    eq=np.diag(f.astype(complex))
    re=_r(e,_c(v,eq),spec.omega_e,spec.eta_e)
    drho=derivative_equilibrium_density(e,u,_data(spec)[1],spec.system.mu,spec.system.temperature)
    deq=np.diag(np.diag(drho))
    off=drho-deq
    # Differentiation of the B-prepared linear E susceptibility:
    # R_E [V,Df(H)[Z]] + R_E [Z, R_E [V,rho0]].
    population=_r(e,_c(v,deq),spec.omega_e,spec.eta_e)
    coherence=_r(e,_c(v,off),spec.omega_e,spec.eta_e)
    e_then_b=_r(e,_c(z,re),spec.omega_e,spec.eta_e)
    val=lambda x:complex(np.trace(a@x))
    return dict(population=val(population),coherence=val(coherence),e_then_b=val(e_then_b),
                total=val(population+coherence+e_then_b),
                static_without_population=val(coherence+e_then_b))


def time_domain_rk4(spec:TwoColorSpec,t_start=-85.,dt=.035):
    """Independent initial-value integration of EACH perturbative Liouville order.

    H0 diagonalized only to speed commutator. Initial rho_E,B,EB(t_start)=0.
    In the remote-past limit these converge to harmonic analytic kernels.
    """
    if t_start>=0 or dt<=0:raise ValueError('invalid time integration range')
    e,z,v,a,f,u=_basis(spec)
    n=len(e);d=e[:,None]-e[None,:];eq=np.diag(f)
    se=_c(v,eq);sb=_c(z,eq)
    steps=int(np.ceil(-t_start/dt));dt=-t_start/steps
    st=np.zeros((3,n,n),dtype=complex)
    s_e=spec.eta_e-1j*spec.omega_e
    s_b=spec.eta_b-1j*spec.omega_b
    def deriv(t,state):
        re,rb,rbe=state
        ee=np.exp(s_e*t);bb=np.exp(s_b*t)
        out=np.empty_like(state)
        out[0]=-1j*(d*re+ee*se)
        out[1]=-1j*(d*rb+bb*sb)
        out[2]=-1j*(d*rbe+ee*_c(v,rb)+bb*_c(z,re))
        return out
    t=t_start
    for _ in range(steps):
        k1=deriv(t,st)
        k2=deriv(t+.5*dt,st+.5*dt*k1)
        k3=deriv(t+.5*dt,st+.5*dt*k2)
        k4=deriv(t+dt,st+dt*k3)
        st=st+(dt/6.)*(k1+2*k2+2*k3+k4)
        t+=dt
    return {'total':complex(np.trace(a@st[2])),
            'rho_e':st[0], 'rho_b':st[1], 'rho_eb':st[2],
            'steps':steps, 't_start':t_start}


def static_eta_scan(spec:TwoColorSpec,etas):
    """Diagnostic, not a tau extrapolation: switches B from adiabatic to prepared."""
    from dataclasses import replace
    if spec.omega_b != 0:raise ValueError('omega_b must be zero')
    prep=static_preparation(spec)
    rows=[]
    for eta in etas:
        p=replace(spec,eta_b=float(eta))
        d=frequency_channels(p)
        rows.append(dict(eta_b=float(eta),dynamic_total=d['total'],
                         difference_from_prepared=d['total']-prep['total'],
                         distance_to_prepared_without_population=abs(d['total']-prep['static_without_population'])))
    return prep,rows


def finite_amplitude_time_domain(spec:TwoColorSpec,amplitude_e=.008,amplitude_b=.006,
                                 t_start=-155.,dt=.045):
    """Independent nonlinear IVP solved at all four signs of E and B.

    The mixed coefficient is the central 4-corner divided difference, so
    unperturbed, E-only, B-only, E^2 and B^2 contributions cancel.
    The remaining O(E^2,B^2) bias is explicitly tested by amplitude halving.
    """
    if amplitude_e<=0 or amplitude_b<=0 or t_start>=0 or dt<=0:
        raise ValueError('positive field amplitudes, t_start<0 and dt>0 required')
    e,z,v,a,f,u=_basis(spec)
    d=e[:,None]-e[None,:]
    eq=np.diag(f.astype(complex))
    se=spec.eta_e-1j*spec.omega_e
    sb=spec.eta_b-1j*spec.omega_b
    nsteps=int(np.ceil(-t_start/dt));step=-t_start/nsteps
    def evolve(ee,bb):
        state=eq.copy();t=t_start
        def deriv(t,mat):
            ve=ee*np.exp(se*t);zb=bb*np.exp(sb*t)
            return -1j*(d*mat+ve*_c(v,mat)+zb*_c(z,mat))
        for _ in range(nsteps):
            k1=deriv(t,state)
            k2=deriv(t+step*.5,state+step*.5*k1)
            k3=deriv(t+step*.5,state+step*.5*k2)
            k4=deriv(t+step,state+step*k3)
            state=state+step*(k1+2*k2+2*k3+k4)/6
            t+=step
        return complex(np.trace(a@state))
    pp=evolve(amplitude_e,amplitude_b)
    pm=evolve(amplitude_e,-amplitude_b)
    mp=evolve(-amplitude_e,amplitude_b)
    mm=evolve(-amplitude_e,-amplitude_b)
    return (pp-pm-mp+mm)/(4*amplitude_e*amplitude_b)
