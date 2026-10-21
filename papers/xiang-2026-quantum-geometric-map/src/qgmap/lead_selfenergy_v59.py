"""v59: Exact surface self energy of semi-infinite 1D TB electrodes.

Finite-device NEGF with energy-dependent (and optionally physical B-dependent)
local contact Σ^r, not invented constant Γ. This is a coherent FINITE DEVICE.
A 1D lead is an exactly specified simple reservoir, NOT a general material lead,
and the resulting Hall voltage is NOT intrinsic bulk orbital MNHE.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import quad_vec
from .embedding_magnetic_v55 import peierls_field
from .lead_transport_v57 import four_terminal, diamond_graph


@dataclass(frozen=True)
class ChainLead:
    site: int
    hopping: float = 2.0
    coupling: float = 0.7416198487095663  # sqrt(0.55) at t=2, Γ(E=0)=0.55
    onsite: float = 0.0
    onsite_B_slope: float = 0.0

    def __post_init__(self):
        if not isinstance(self.site, (int, np.integer)) or self.site < 0:
            raise ValueError('site must be a nonnegative integer')
        if not np.isfinite([self.hopping,self.coupling,self.onsite,self.onsite_B_slope]).all() or self.hopping<=0 or self.coupling<=0:
            raise ValueError('physical semi-infinite lead needs finite t>0 and |w|>0')


def surface_g(E, t, onsite=0):
    """Retarded g_11(E+i0) of a *semi-infinite* uniform 1D nearest-neighbour lead.

    Branch is analytic upper half-plane, g~1/E; spectral Γ >=0. Works also
    outside the lead band, where Im g=0, and exactly at band edges.
    """
    if not np.isfinite([E,t,onsite]).all() or t<=0:raise ValueError('real finite E, onsite and t>0')
    x=E-onsite
    edge=2*t
    if x>=edge:
        return complex(2/(x+np.sqrt(x*x-edge*edge)))
    if x<=-edge:
        return complex(2/(x-np.sqrt(x*x-edge*edge)))
    return complex((x-1j*np.sqrt(max(0.,edge*edge-x*x)))/(2*t*t))


def surface_g_dx(E,t,onsite=0):
    """Analytic energy derivative from g=(E-onsite-t²g)^-1.

    d g / dE= -g²/(1-t²g²); divergent exactly at band edges. Our energy
    quadrature never evaluates endpoints as interior nodes. Flag edges.
    """
    g=surface_g(E,t,onsite)
    den=1-t*t*g*g
    if abs(den)<1e-11:raise ValueError('surface derivative divergent at 1D band edge')
    return -g*g/den


def selfenergy_and_B(E,B,lead:ChainLead):
    """Explicit magnetic lead effect is optional: onsite(B)=onsite+slope*B.

    Σ=w²g(E-onsite(B)), d_BΣ=-slope*w² dg/dE.
    Independent physical input. NOT a Peierls magnetic gauge transform.
    """
    center=lead.onsite+lead.onsite_B_slope*B
    sigma=lead.coupling**2*surface_g(E,lead.hopping,center)
    if lead.onsite_B_slope==0:
        deriv=0j
    else:
        deriv=-lead.onsite_B_slope*lead.coupling**2*surface_g_dx(E,lead.hopping,center)
    gamma=-2*sigma.imag
    dgamma=-2*deriv.imag
    if gamma< -1e-12:raise RuntimeError('causality violation: negative Gamma')
    return sigma,deriv,max(0.,gamma),dgamma


@dataclass(frozen=True)
class ChainSpec:
    field:float=0.31
    charge_over_hbar:float=-.87
    mu:float=.15
    temperature:float=.16
    thermal_window:float=18.
    epsrel:float=1e-8
    epsabs:float=1e-10
    def __post_init__(self):
        v=(self.field,self.charge_over_hbar,self.mu,self.temperature,self.thermal_window,self.epsrel,self.epsabs)
        if not np.isfinite(v).all() or self.charge_over_hbar==0 or self.temperature<=0 or self.thermal_window<10 or self.epsrel<=0 or self.epsabs<=0:
            raise ValueError('finite physical E/B, positive T/integration tolerances')


def energy_transmissions(E,H,M,leads,B):
    """Exact Caroli T_ab(E,B) and total d_B T incl. Σ'(E,B) and Γ'(E,B).

    Valid for B-independent device-lead hopping and diagonal local reservoir
    projectors. When lead onsite shifts with B, both Σ' and Γ' are retained.
    """
    H=np.asarray(H,complex);M=np.asarray(M,complex)
    n=H.shape[0]
    if H.shape!=(n,n) or M.shape!=(n,n) or n<2 or len(leads) not in (2,4):raise ValueError('matrix/contact shapes')
    if not np.allclose(H,H.conj().T,atol=1e-12) or not np.allclose(M,M.conj().T,atol=1e-12):raise ValueError('Hermitian H, M required')
    sites=[lead.site for lead in leads]
    if len(set(sites))!=len(sites) or any(s>=n for s in sites):raise ValueError('distinct valid terminal sites')
    nleads=len(leads)
    selfener=np.zeros(n,complex);dself=np.zeros(n,complex)
    gammas=np.zeros(nleads);dgammas=np.zeros(nleads)
    for i,lead in enumerate(leads):
        selfener[lead.site],dself[lead.site],gammas[i],dgammas[i]=selfenergy_and_B(E,B,lead)
    # With all physically open channels closed, transmission is identically
    # zero; avoid spurious singular matrix at a decoupled bound-state energy.
    if np.sum(gammas)<1e-14:
        return np.zeros((nleads,nleads)),np.zeros((nleads,nleads))
    gr=np.linalg.inv(E*np.eye(n,dtype=complex)-H-np.diag(selfener))
    dgr=gr@(M+np.diag(dself))@gr
    T=np.zeros((nleads,nleads));dT=T.copy()
    for a in range(nleads):
        for b in range(nleads):
            if a==b:continue
            g=gr[sites[a],sites[b]]
            dg=dgr[sites[a],sites[b]]
            strength=gammas[a]*gammas[b]
            dstrength=dgammas[a]*gammas[b]+gammas[a]*dgammas[b]
            T[a,b]=strength*abs(g)**2
            dT[a,b]=dstrength*abs(g)**2+2*strength*np.real(g.conjugate()*dg)
    return T,dT


def integrate_leads(H0,positions,leads,spec=ChainSpec(),*,gauge='symmetric'):
    """Adaptive full-band Fermi-window Caroli integration, with band edges.

    Finite-temperature linear-response weights, with explicit integration
    error, optical-theorem gate, and all floating voltage contact derivatives.
    """
    h0=np.asarray(H0,complex);positions=np.asarray(positions,float)
    if h0.ndim!=2 or h0.shape[0]!=h0.shape[1] or positions.shape!=(len(h0),2):raise ValueError('graph arrays')
    pb=peierls_field(h0,positions,spec.field,spec.charge_over_hbar,gauge=gauge)
    low=spec.mu-spec.thermal_window*spec.temperature
    high=spec.mu+spec.thermal_window*spec.temperature
    # Independently derived 1D lead band edges; also split at finite eigenvalues
    # to capture narrow mesoscopic resonances.
    edges=[]
    for l in leads:
        c=l.onsite+l.onsite_B_slope*spec.field
        edges.extend([c-2*l.hopping,c+2*l.hopping])
    edges.extend(np.linalg.eigvalsh(pb.hamiltonian))
    points=sorted(set(float(e) for e in edges if low<e<high))
    count=[0]
    def integrand(E):
        count[0]+=1
        tt,dt=energy_transmissions(E,pb.hamiltonian,pb.magnetic_vertex,leads,spec.field)
        f=1/(4*spec.temperature*np.cosh((E-spec.mu)/(2*spec.temperature))**2)
        return np.r_[tt.ravel(),dt.ravel()]*f
    m=len(leads)
    vec,error,info=quad_vec(integrand,low,high,points=points,epsrel=spec.epsrel,epsabs=spec.epsabs,limit=700,full_output=True)
    if not info.success:raise RuntimeError('nonconverged energy-integrated lead response')
    T=vec[:m*m].reshape(m,m);dT=vec[m*m:].reshape(m,m)
    # Near-local real lead Hamiltonian ensures exact conservation in the
    # complete energy integral (no absorption / phenomenological baths).
    curr_residual=float(np.max(abs(T.sum(axis=0)-T.sum(axis=1))))
    dcurr_residual=float(np.max(abs(dT.sum(axis=0)-dT.sum(axis=1))))
    if curr_residual>max(2e-7,30*error) or dcurr_residual>max(2e-6,100*error):
        raise RuntimeError('lead current conservation broken; check spectral weights and contacts')
    result={'T':T,'dT_dB':dT,'quadrature_error':float(error),'quadrature_evaluations':count[0],
            'current_conservation_residual':curr_residual,'d_current_conservation_residual':dcurr_residual,
            'status':'VERIFIED_SPECIFIED_SEMIINFINITE_1D_LEAD_ALGORITHM_CONDITIONAL_DEVICE'}
    if m==4:result['observables']=four_terminal(T,dT)
    else:result['observables']={'two_terminal_conductance_q2_over_h':float(T[0,1]),
                                'dTwo_terminal_dB':float(dT[0,1])}
    return result


def finite_B_fd(H0,positions,leads,spec=ChainSpec(),*,gauge='symmetric',step=2e-4):
    """Independent 5-point full device differentiation; includes Σ(B) and Γ(B)."""
    if not np.isfinite(step) or step<=0:raise ValueError('positive B step')
    vals=[integrate_leads(H0,positions,leads,replace(spec,field=spec.field+d*step),gauge=gauge) for d in (2,1,-1,-2)]
    def five(get):return (-get(vals[0])+8*get(vals[1])-8*get(vals[2])+get(vals[3]))/(12*step)
    out={'dT_dB':five(lambda z:z['T'])}
    if len(leads)==4:
        out['dHall_dB']=float(five(lambda z:z['observables']['hall_ratio']))
    else:out['dTwo_terminal_dB']=float(five(lambda z:z['observables']['two_terminal_conductance_q2_over_h']))
    return out


def reference_device(spec=ChainSpec(), lead_hopping=2., onsite_B_slope=0.):
    pos,H0,sites=diamond_graph(1)
    # Choose w to make Γ(E=onsite)=0.55 at B=0, so wide-band t->infty
    # limit has fixed 0.55, not increasing lead strength.
    w=np.sqrt(.55*lead_hopping/2)
    leads=tuple(ChainLead(site=s,hopping=lead_hopping,coupling=w,onsite_B_slope=onsite_B_slope) for s in sites)
    return H0,pos,leads
