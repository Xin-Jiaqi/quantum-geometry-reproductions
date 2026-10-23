"""v75 explicitly elastic, charge-conserving Dirac angular collision gate.

Independent physical precursor for the intraband electric lifetime. The
isothermal spin Zeeman E*B current of Xiang PRB is NOT generated here.
The scalar impurity overlap (1+cos dtheta)/2 follows massless Dirac
helicity. The rate gamma is an *input* Born-style angular scattering scale,
not a rate extracted from any Wannier or experimental disorder data.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Spec:
    beta:float=.125
    energy:float=.3
    gamma:float=.07
    n:int=96
    spin_overlap:float=1.
    def __post_init__(self):
        vals=(self.beta,self.energy,self.gamma,self.spin_overlap)
        if not np.isfinite(vals).all() or abs(self.beta)>=1 or self.energy<=0 or self.gamma<=0 or not 0<=self.spin_overlap<=1:
            raise ValueError('type-I Dirac tilt, energy/gamma>0, 0<=overlap<=1 required')
        if not isinstance(self.n,int) or self.n<12 or self.n%4:
            raise ValueError('n multiple of four and >=12 required')


def collision_data(s=Spec()):
    theta=2*np.pi*np.arange(s.n)/s.n
    cost=np.cos(theta);sint=np.sin(theta)
    denom=1+s.beta*cost
    r=s.energy/denom
    eps=r*denom
    weight=denom**-2
    weight=weight/np.sum(weight)
    # Scalar impurity spinor overlap: 2|<u(theta)|u(theta')>|^2=1+cos(theta-theta').
    # With detailed-balance final-state density weights w_j.
    matrix=s.gamma*weight[None,:]*(1+s.spin_overlap*np.cos(theta[:,None]-theta[None,:]))
    L=np.diag(matrix.sum(axis=1))-matrix
    # For g = delta f / (-f'(energy)), L g = E v_y.
    vy=sint.copy()
    vx=s.beta+cost
    return {'theta':theta,'radius':r,'energies':eps,'weights':weight,'collision':matrix,
            'L':L,'vy':vy,'vx':vx,'row_rates':matrix.sum(axis=1)}


def steady(s=Spec(),omega=0.):
    if not np.isfinite(omega):raise ValueError('finite frequency required')
    d=collision_data(s)
    w=d['weights'];v=d['vy'];L=d['L']
    v0=v-np.dot(w,v)
    # Fix the conserved zero-mode by a number-conserving projector. No
    # fictitious density-decay gamma or full-sector broadening is introduced.
    sol=np.linalg.solve(L-1j*omega*np.eye(s.n)+np.outer(np.ones(s.n),w),v0)
    residual=(L-1j*omega*np.eye(s.n))@sol-v0
    curr=np.dot(w,v*sol)
    tau=curr/np.dot(w,v*v)
    return {'distribution':sol,'tau_current':tau,'current':curr,
            'weighted_number':np.dot(w,sol),'collision_charge_residual':np.max(abs(w@L)),
            'linear_solver_residual':np.max(abs(residual)),'source_weighted_mean':np.dot(w,v0),
            'drude_amplitude':np.dot(w,v*v),'status':'ELASTIC_NUMBER_CONSERVING_SINGLE_BAND_NOT_FULL_PHE'}


def eig_physics(s=Spec()):
    d=collision_data(s);w=d['weights'];L=d['L']
    D=np.sqrt(w)
    sym=D[:,None]*L/D[None,:]
    z=np.linalg.eigvalsh((sym+sym.T)*.5)
    return {'eigenvalues':z,'zero_mode':z[0],'spectral_gap':z[1],
            'detailed_balance_error':np.max(abs(w[:,None]*d['collision']-w[None,:]*d['collision'].T)),
            'energy_shell_error':np.max(abs(d['energies']-s.energy)),
            'number_conservation_error':np.max(abs(w@L)),
            'constant_mode_error':np.max(abs(L@np.ones(s.n))),
            'equal_radius_energy_spread':float(2*abs(s.beta)*s.energy),
            'status':'BORN_STYLE_ELASTIC_ISOENERGY_KERNEL_NOT_DISORDER_FIT'}


def sign_invariance(energies,position,spin,df):
    """Source Eq20 -> Eq31 *relative sign* cannot change by common conventions.

    This applies to source v49/v71 literal epsilon_nm and velocity identity,
    not to a hypothetical alternative printed definition of Z or spin source.
    """
    from .source_observable_bridge_v71 import source_eq20_eq31_sign
    ref=source_eq20_eq31_sign(energies,position,spin,df)
    flip_s=source_eq20_eq31_sign(energies,position,-spin,df)
    flip_r=source_eq20_eq31_sign(energies,-position,spin,df)
    return {'literal':ref,'flipped_spin':flip_s,'flipped_position':flip_r,
            'relative_sum':ref['direct_eq20']+ref['printed_eq31'],
            'spin_flip_both':(flip_s['direct_eq20']+ref['direct_eq20'],
                              flip_s['printed_eq31']+ref['printed_eq31']),
            'position_flip_both':(flip_r['direct_eq20']+ref['direct_eq20'],
                                  flip_r['printed_eq31']+ref['printed_eq31']),
            'status':'LITERAL_SOURCE_RELATIVE_SIGN_IRREDUCIBLE_WITH_COMMON_CONVENTION_FLIP'}


def independent_weighted_spectral(s=Spec(),omega=0.):
    """Independent weighted-Hermitian diagonalization of detailed-balance L.

    Removes the conserved charge mode; direct matrix inversion is NOT reused.
    """
    if not np.isfinite(omega):raise ValueError('finite omega required')
    d=collision_data(s)
    w=d['weights'];wroot=np.sqrt(w);v=d['vy']
    S=wroot[:,None]*d['L']/wroot[None,:]
    ev,Q=np.linalg.eigh((S+S.T)/2)
    v0=v-np.dot(w,v)
    coeff=Q.T@(wroot*v0)
    coeff[np.abs(ev)<1e-10]=0
    p=np.divide(coeff,ev-1j*omega,out=np.zeros(len(ev),complex),where=np.abs(ev)>1e-10)
    g=(Q@p)/wroot
    return {'distribution':g,'tau_current':np.dot(w,v*g)/np.dot(w,v*v),
            'zero_mode_count':int(np.count_nonzero(abs(ev)<1e-10)),
            'positive_eigenvalue_min':float(ev[1])}


def time_domain_independent(s=Spec(n=24),duration=200.):
    """Independent driven finite-time solution, with constant source and g(0)=0.

    dg/dt = v - Lg, and after t the response is integral exp(-Ls) v ds.
    """
    from scipy.linalg import expm
    if not np.isfinite(duration) or duration<=0:raise ValueError('positive duration required')
    d=collision_data(s);w=d['weights'];L=d['L']
    v=d['vy']-np.dot(w,d['vy'])
    # Add projector to invert only the zero-mode-free space.
    P=np.outer(np.ones(s.n),w)
    finite=np.linalg.solve(L+P,(np.eye(s.n)-expm(-L*duration))@v)
    return {'distribution':finite,'weighted_number':float(np.real(np.dot(w,finite)))}
