"""v73: two-color (scalar-potential E_q, uniform spin-Zeeman B) Bloch response.

Retains BOTH E->B, B->E orderings and separate retarded complex poles
w_e=omega_e+i eta_e, w_b=omega_b+i eta_b, w_sum=w_e+w_b.
Finite-q scalar potential is *not* a homogeneous length gauge. A prepared
static B population channel is separated from the collisionless two-color
response; it is not an author's tau-pole or a source-accurate steady state.
"""
from dataclasses import dataclass,replace
import numpy as np
from .bulk_finite_q_spin_v72 import (
    Spec as ParentSpec, _entries, _dag, _tr, fermi, _frechet_f,
    _response, analytic_zeeman_mixed)

@dataclass(frozen=True)
class Spec:
    n:int=24
    q_steps:int=1
    omega_e:float=.27
    omega_b:float=.12
    eta_e:float=.045
    eta_b:float=.03
    B0:float=.11
    mu:float=.32
    temperature:float=.12
    beta:float=.24
    alpha:float=.13
    mass:float=.32
    wilson:float=.35
    orbital_gauge:float=0.
    def __post_init__(self):
        vals=(self.omega_e,self.omega_b,self.eta_e,self.eta_b,self.B0,
              self.mu,self.temperature,self.beta,self.alpha,self.mass,
              self.wilson,self.orbital_gauge)
        if not np.isfinite(vals).all() or self.eta_e<=0 or self.eta_b<=0 or self.temperature<=0:
            raise ValueError('positive finite eta_e, eta_b, temperature required')
        if not isinstance(self.n,int) or self.n<8 or not isinstance(self.q_steps,int) or not 1<=self.q_steps<self.n//2:
            raise ValueError('commensurate nonzero finite q required')

def _parent(s):
    return ParentSpec(n=s.n,q_steps=s.q_steps,omega=s.omega_e,eta=s.eta_e,
        field=s.B0,mu=s.mu,temperature=s.temperature,beta=s.beta,
        alpha=s.alpha,mass=s.mass,wilson=s.wilson,
        orbital_gauge=s.orbital_gauge)

def _eig_data(s):
    hp,hm,zp,zm,W,Jx,Jy,q=_entries(_parent(s))
    ep,up=np.linalg.eigh(hp);em,um=np.linalg.eigh(hm)
    fp=fermi(ep,s.mu,s.temperature);fm=fermi(em,s.mu,s.temperature)
    Fp=(up*fp[...,None,:])@_dag(up)
    Fm=(um*fm[...,None,:])@_dag(um)
    return locals()

def _spectral_resolvent(ep,em,up,um,source,w):
    den=w-(ep[..., :,None]-em[...,None,:])
    out=(_dag(up)@source@um)/den
    return up@out@_dag(um)

def _population_matrix(e,u,z,s):
    # Diagonal equal-energy Fermi derivative f'(epsilon)*Z_nn.
    occ=fermi(e,s.mu,s.temperature)
    zband=_dag(u)@z@u
    diagonal=-occ*(1-occ)/s.temperature*np.diagonal(zband,axis1=-2,axis2=-1)
    bd=np.zeros_like(zband)
    idx=np.arange(2)
    bd[...,idx,idx]=diagonal
    return u@bd@_dag(u)

def _source_and_currents(R,D,w):
    density,lon,tra=_response(R,D['W'],D['Jx'],D['Jy'])
    return dict(density=density,longitudinal=lon,transverse=tra,
                normalized_E_transverse=tra/(-1j*D['q']))

def two_color(s=Spec(),prepared_static_population=False):
    """Both time orderings; optional static equilibrium population preparation.

    Collisionless source: RB=L_B^-1([Z,F]), RE=L_E^-1(WF- -F+W).
    BE=L_sum^-1(W RB- - RB+ W), EB=L_sum^-1(Z+RE-RE Z-).
    For omega_b=0, only by adding a separately declared f'(H)Zdiag
    channel to RB can eta_b->0 recover the fully equilibrated d_B RE.
    This channel is a preparation prescription, not source SM Eq15's pole.
    """
    if prepared_static_population and s.omega_b!=0.:
        raise ValueError('equilibrium population preparation requires omega_b=0')
    D=_eig_data(s)
    hp,hm,zp,zm,W=D['hp'],D['hm'],D['zp'],D['zm'],D['W']
    ep,em,up,um=D['ep'],D['em'],D['up'],D['um']
    Fp,Fm=D['Fp'],D['Fm']
    we=s.omega_e+1j*s.eta_e
    wb=s.omega_b+1j*s.eta_b
    ws=we+wb
    source_e=W@Fm-Fp@W
    RE=_spectral_resolvent(ep,em,up,um,source_e,we)
    RBp=_spectral_resolvent(ep,ep,up,up,zp@Fp-Fp@zp,wb)
    RBm=_spectral_resolvent(em,em,um,um,zm@Fm-Fm@zm,wb)
    if prepared_static_population:
        RBp=RBp+_population_matrix(ep,up,zp,s)
        RBm=RBm+_population_matrix(em,um,zm,s)
    source_BE=W@RBm-RBp@W
    source_EB=zp@RE-RE@zm
    RBE=_spectral_resolvent(ep,em,up,um,source_BE,ws)
    REB=_spectral_resolvent(ep,em,up,um,source_EB,ws)
    Rtotal=RBE+REB
    observables={k:_source_and_currents(R,D,ws) for k,R in
                 [('b_then_e',RBE),('e_then_b',REB),('total',Rtotal)]}
    rhs=complex(np.mean(_tr(_dag(W)@(source_BE+source_EB))))
    density,lon=observables['total']['density'],observables['total']['longitudinal']
    ward=ws*density-D['q']*lon-rhs
    return dict(**observables,q=D['q'],total_source_integrated=rhs,
                ward_residual=complex(ward),frequency_sum=ws,
                electric_pole=we,magnetic_pole=wb,
                prepared_static_population=bool(prepared_static_population),
                status=('PREPARED_STATIC_POPULATION_SEPARATE_PROTOCOL_NOT_SOURCE_TAU_POLE'
                        if prepared_static_population else
                        'TWO_COLOR_FINITE_Q_COLLISIONLESS_SPIN_NOT_SOURCE_PHE_DC'))

def _dense_resolvent(A,B,S,w):
    """Independent 2x2 orbital Kronecker Sylvester inversion at each Bloch k."""
    from .bulk_finite_q_spin_v72 import I
    out=np.empty_like(S)
    for a in range(len(S)):
        for b in range(len(S[a])):
            linear=np.kron(w*I-A[a,b],I)+np.kron(I,B[a,b].T)
            out[a,b]=np.linalg.solve(linear,S[a,b].reshape(4)).reshape(2,2)
    return out

def two_color_dense_independent(s=Spec(),prepared_static_population=False):
    if prepared_static_population and s.omega_b!=0.:raise ValueError('static only')
    D=_eig_data(s)
    hp,hm,zp,zm,W=D['hp'],D['hm'],D['zp'],D['zm'],D['W']
    Fp,Fm=D['Fp'],D['Fm']
    we=s.omega_e+1j*s.eta_e
    wb=s.omega_b+1j*s.eta_b
    ws=we+wb
    RE=_dense_resolvent(hp,hm,W@Fm-Fp@W,we)
    RBp=_dense_resolvent(hp,hp,zp@Fp-Fp@zp,wb)
    RBm=_dense_resolvent(hm,hm,zm@Fm-Fm@zm,wb)
    if prepared_static_population:
        RBp+=_population_matrix(D['ep'],D['up'],zp,s)
        RBm+=_population_matrix(D['em'],D['um'],zm,s)
    RBE=_dense_resolvent(hp,hm,W@RBm-RBp@W,ws)
    REB=_dense_resolvent(hp,hm,zp@RE-RE@zm,ws)
    return {k:_source_and_currents(R,D,ws) for k,R in
            [('b_then_e',RBE),('e_then_b',REB),('total',RBE+REB)]}

def static_magnetic_derivative(s=Spec(omega_b=0.,eta_b=1e-5)):
    """Compare collisionless vs prepared-two-color static B to exact dB response."""
    if s.omega_b!=0.:raise ValueError('requires static magnetic frequency')
    ref=analytic_zeeman_mixed(_parent(s))
    collisionless=two_color(s,False)
    prepared=two_color(s,True)
    name='normalized_E_transverse'
    target=ref['dB_E_transverse']
    return dict(exact_prepared_dB=target,
                collisionless=collisionless['total'][name],
                with_population=prepared['total'][name],
                collisionless_error=float(abs(collisionless['total'][name]-target)),
                prepared_error=float(abs(prepared['total'][name]-target)),
                missing_population_signal=complex(prepared['total'][name]-collisionless['total'][name]),
                source_limit_status='FINITE_Q_PREPARED_B_LIMIT_NOT_SELECTED_TAU_POLE')

def regulator_order_scan(s=Spec()):
    """Numerical diagnosis only: finite-q pole and dc uniform E not identified."""
    if s.omega_e!=0.:
        s=replace(s,omega_e=0.)
    vals={}
    for eta in (.1,.05,.025,.0125):
        ss=replace(s,eta_e=eta,eta_b=.75*eta)
        v=two_color(ss)['total']['normalized_E_transverse']
        vals[str(eta)]={'value':v,'eta_e_times_value':eta*v}
    return dict(values=vals,status='NO_SOURCE_TAU_IDENTIFICATION_AT_FIXED_FINITE_Q')
