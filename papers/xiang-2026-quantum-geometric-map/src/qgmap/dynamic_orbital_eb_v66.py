"""v66 finite-frequency projected-orbital EB with ALL magnetic + electric contacts.

A finite OPEN Peierls graph and SUPPLIED projected P_x/P_y matrices. This is a
conditional scalar-gamma Liouville bath, NOT conserving Kubo--Streda, not the
Xiang source's intrinsic dc orbital magnetotransport.

H_E=H_B-E_b P_b(B), J_a(E)=i[H_E,P_a(B)]. Projected P_a,P_b need not commute,
therefore sigma_ab includes the electric current contact -i[P_b,P_a].
B derivative includes all: d_B J, d_B P_b in drive, rho'_eq, d_B
Liouvillian, d_B electric contact times rho0, and contact times rho'_eq.
"""
from dataclasses import dataclass
from dataclasses import replace
import numpy as np
from .equilibrium_orbital_eb_v65 import (
    Spec as StaticSpec, operators, finite_fermi_frechet,fermi, benchmark_graph,
    full_projected_polarization)


def comm(a,b): return a@b-b@a


def _resolvent(H,source,z):
    if z.real<=0 or not np.isfinite([z.real,z.imag]).all(): raise ValueError('positive gamma required')
    e,U=np.linalg.eigh(H);su=U.conj().T@source@U
    return U@(su/(z+1j*(e[:,None]-e[None,:])))@U.conj().T


@dataclass(frozen=True)
class Spec:
    mu:float=.13
    temperature:float=.37
    field:float=.23
    charge_over_hbar:float=-.84
    gamma:float=.071
    omega:float=.29
    def __post_init__(self):
        if not np.isfinite([self.mu,self.temperature,self.field,self.charge_over_hbar,self.gamma,self.omega]).all():
            raise ValueError('finite physical parameters')
        if self.temperature<=0 or self.gamma<=0 or self.charge_over_hbar==0:
            raise ValueError('positive temperature,gamma; nonzero charge')


def _static_spec(s):return StaticSpec(s.mu,s.temperature,s.field,s.charge_over_hbar)


def _ops(hoppings,positions,Pa,Pb,spec,frame,K):
    sa=operators(hoppings,positions,Pa,_static_spec(spec),frame,K)
    sb=operators(hoppings,positions,Pb,_static_spec(spec),frame,K)
    if not np.allclose(sa['H'],sb['H']) or not np.allclose(sa['M'],sb['M']):
        raise RuntimeError('field/operator mismatch')
    return sa['H'],sa['M'],sa['P'],sa['dP'],sb['P'],sb['dP']


def ac_state(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None, electric_contact=True):
    """Finite-omega E response with optional omitted-contact negative control."""
    H,M,A,dA,B,dB=_ops(hoppings,positions,Pa,Pb,spec,frame,K)
    eig,U=np.linalg.eigh(H)
    rho=(U*fermi(eig,spec.mu,spec.temperature)[None,:])@U.conj().T
    J=1j*comm(H,A)
    JE=-1j*comm(B,A)
    src=1j*comm(B,rho)
    z=spec.gamma-1j*spec.omega
    rhoE=_resolvent(H,src,z)
    sig=complex(np.trace(J@rhoE)+(np.trace(JE@rho) if electric_contact else 0))
    return {'sigma':sig,'H':H,'M':M,'Pa':A,'dPa':dA,'Pb':B,'dPb':dB,
            'rho':rho,'rhoE':rhoE,'J':J,'JE':JE,'source':src,'z':z,
            'status':'FINITE_PROJECTED_SCALAR_GAMMA_AC_CONDITIONAL_NOT_ORIGINAL_ORBITAL_EB'}


def dynamic_mixed_eb(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None, electric_contact=True):
    """Exact analytic B derivative of complete AC response; seven named parts.

    Source d_B has separate drive-position and equilibrium preparation terms;
    JE d_B has both dPa/dPb; J d_B has H and Pa derivatives.
    """
    s=ac_state(hoppings,positions,Pa,Pb,spec,frame,K,electric_contact)
    H,M,A,dA,B,dB=s['H'],s['M'],s['Pa'],s['dPa'],s['Pb'],s['dPb']
    rho,rhoE,J,JE=s['rho'],s['rhoE'],s['J'],s['JE']
    z=s['z']
    rhoB=finite_fermi_frechet(H,M,spec.mu,spec.temperature)
    J_H=1j*comm(M,A)
    J_P=1j*comm(H,dA)
    source_P=_resolvent(H,1j*comm(dB,rho),z)
    source_rho=_resolvent(H,1j*comm(B,rhoB),z)
    dynamics=_resolvent(H,-1j*comm(M,rhoE),z)
    JE_B=-1j*(comm(dB,A)+comm(B,dA))
    components={
        'current_magnetic_H':complex(np.trace(J_H@rhoE)),
        'current_position_dPa':complex(np.trace(J_P@rhoE)),
        'source_position_dPb':complex(np.trace(J@source_P)),
        'source_fermi_density':complex(np.trace(J@source_rho)),
        'liouvillian_magnetic':complex(np.trace(J@dynamics)),
        'electric_current_contact_dB':complex(np.trace(JE_B@rho)) if electric_contact else 0j,
        'electric_current_contact_rhoB':complex(np.trace(JE@rhoB)) if electric_contact else 0j,
    }
    return {'sigma':s['sigma'],'mixed_dsigma_dB':sum(components.values()),
            'components':components,'frame':frame,
            'equilibrium_dB_density_trace':complex(np.trace(rhoB)),
            'status':'CONDITIONAL_FINITE_AC_MIXED_EB_GAMMA_BATH_NOT_CONSERVING_STREDA'}


def ac_spectral_independent(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None, electric_contact=True):
    """Independent explicit double-eigenstate Kubo-like sum, no resolvent method."""
    H,M,A,dA,B,dB=_ops(hoppings,positions,Pa,Pb,spec,frame,K)
    e,U=np.linalg.eigh(H)
    f=fermi(e,spec.mu,spec.temperature)
    current=U.conj().T@(1j*comm(H,A))@U
    dipole=U.conj().T@B@U
    value=0j
    for m in range(len(e)):
        for n in range(len(e)):
            value+=current[n,m]*(1j*dipole[m,n]*(f[n]-f[m]))/(spec.gamma-1j*spec.omega+1j*(e[m]-e[n]))
    if electric_contact:
        eq=(U*f[None,:])@U.conj().T
        value+=np.trace((-1j*comm(B,A))@eq)
    return complex(value)


def full_field_finite_difference(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None,
                                 step=1e-4,points=5,electric_contact=True):
    """Completely rediagonalize at five distinct magnetic fields (independent)."""
    if not np.isfinite(step) or step<=0 or points not in (3,5):raise ValueError('positive step; 3/5')
    def response(field):
        return ac_state(hoppings,positions,Pa,Pb,replace(spec,field=field),frame,K,electric_contact)['sigma']
    x=spec.field
    if points==3:return (response(x+step)-response(x-step))/(2*step)
    return (response(x-2*step)-8*response(x-step)+8*response(x+step)-response(x+2*step))/(12*step)


def benchmark_model(tr_breaking=True):
    """Explicitly supplied projected x/y dipoles; both interorbital and Hermitian."""
    h,xy,dy=benchmark_graph(tr_breaking)
    dx=np.zeros_like(dy)
    for i,j,z in ((0,2,.11+.03j),(1,3,-.09+.1j),(0,1,.08-.04j)):
        if not tr_breaking:z=z.real
        dx[i,j]=z;dx[j,i]=z.conjugate()
    s=Spec()
    px=full_projected_polarization(xy,s.charge_over_hbar,'x',dx)
    py=full_projected_polarization(xy,s.charge_over_hbar,'y',dy)
    k=np.array([[.3,.21+.1j,0,0],[.21-.1j,-.4,.11,0],
                [0,.11,.2,.17+.08j],[0,0,.17-.08j,-.1]],complex)
    return h,xy,px,py,k


def independent_dc_steady_current(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None,E=0.):
    """Independent finite-E DC bath steady current via *dense superoperator* solve.

    (gamma+i ad_{H(E)})rho=gamma*f(H_B); current i[H(E),Pa].
    Bath targets field-equilibrated B-only f(H_B), not f(H(B)-E Pb).
    This checks the algebraic omega=0 limit but is NOT microscopic scattering.
    """
    if not np.isfinite(E):raise ValueError('finite E')
    H,_,A,_,B,_=_ops(hoppings,positions,Pa,Pb,spec,frame,K)
    e,U=np.linalg.eigh(H)
    rho0=(U*fermi(e,spec.mu,spec.temperature)[None,:])@U.conj().T
    HE=H-E*B
    n=len(H)
    # column-major vec(H rho-rho H) = (I⊗H-H^T⊗I) vec(rho)
    superop=spec.gamma*np.eye(n*n)+1j*(np.kron(np.eye(n),HE)-np.kron(HE.T,np.eye(n)))
    rho=np.linalg.solve(superop,spec.gamma*rho0.reshape(-1,order='F')).reshape((n,n),order='F')
    J=1j*comm(HE,A)
    return complex(np.trace(J@rho))


def independent_dc_mixed_fd(hoppings,positions,Pa,Pb,spec=Spec(),frame='symmetric',K=None,step=2e-4):
    """Four-corner simultaneous finite B and finite E, no linear-response solver."""
    if not np.isfinite(step) or step<=0:raise ValueError('positive step')
    if abs(spec.omega)>1e-15:raise ValueError('strict DC reference requires omega=0')
    def I(b,e):return independent_dc_steady_current(hoppings,positions,Pa,Pb,replace(spec,field=b),frame,K,e)
    b=spec.field;s=step
    return (I(b+s,s)-I(b+s,-s)-I(b-s,s)+I(b-s,-s))/(4*s*s)
