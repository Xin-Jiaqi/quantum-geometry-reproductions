"""v46: two explicitly different scalar-RTA mixed E/B protocols for tilted Dirac.

Same H(k,b)=beta*kx I+kx sy-(ky+b)sx in both protocols. Source uses
Fermi projector divided derivatives and the full physical-spin Liouvillian.
B-prepared collision target f(H(B)) gives a complete B-derivative constrained
by ky-translation; frozen-H0 collision target yields a different mixed E*B.
This is not an identification with the original paper's supplemental RTA.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss

I=np.eye(2,dtype=complex)
SX=np.array([[0.,1.],[1.,0.]],complex)
SY=np.array([[0.,-1j],[1j,0.]],complex)

@dataclass(frozen=True)
class Spec:
    beta:float=.125
    mu:float=.30
    temperature:float=.055
    eta:float=.065
    cutoff:float=2.5
    nr:int=64
    nphi:int=96
    def __post_init__(self):
        if (not np.isfinite([self.beta,self.mu,self.temperature,self.eta,self.cutoff]).all()
            or abs(self.beta)>=1 or min(self.temperature,self.eta,self.cutoff)<=0
            or self.nr<6 or self.nphi<16):
            raise ValueError('invalid scalar-RTA configuration')

@lru_cache(maxsize=64)
def quadrature(nr,nphi,cutoff):
    q,w=leggauss(nr)
    r=(q+1)*cutoff/2
    wr=w*cutoff/2
    ang=2*np.pi*(np.arange(nphi)+.357)/nphi
    kx=(r[:,None]*np.cos(ang)[None,:]).ravel()
    ky=(r[:,None]*np.sin(ang)[None,:]).ravel()
    wt=np.broadcast_to((wr*r)[:,None]/(2*np.pi*nphi),(nr,nphi)).ravel()
    return kx,ky,wt

def _blocks(spec,kx,ky,field):
    qy=ky+field
    rr=np.hypot(kx,qy)
    if np.any(rr<1e-13):raise ValueError('Dirac degeneracy at quadrature point')
    es0=spec.beta*kx
    en=np.stack((es0+rr,es0-rr),axis=1)
    nmat=(-qy/rr)[:,None,None]*SX+(kx/rr)[:,None,None]*SY
    pp=(I+nmat)*.5;pm=(I-nmat)*.5
    f=1/(1+np.exp(np.clip((en-spec.mu)/spec.temperature,-500,500)))
    fp=-f*(1-f)/spec.temperature
    df=np.stack((qy/rr,-qy/rr),axis=1)
    # v47 correction: d(ky/r)/dky has opposite sign to d(kx/r)/dkx.
    # Independently gated against centered finite differences of f(H_B).
    dn=(-kx*kx/rr**3)[:,None,None]*SX-(kx*qy/rr**3)[:,None,None]*SY
    dpp=dn*.5
    rhoeq=f[:,0,None,None]*pp+f[:,1,None,None]*pm
    source=(fp[:,0]*df[:,0])[:,None,None]*pp+(fp[:,1]*df[:,1])[:,None,None]*pm+(f[:,0]-f[:,1])[:,None,None]*dpp
    return en,pp,pm,rhoeq,source

def _inverse(spec,en,pp,pm,source):
    gap=en[:,0]-en[:,1];eta=spec.eta
    return ((pp@source@pp+pm@source@pm)/eta
      +(pp@source@pm)/(eta+1j*gap)[:,None,None]
      +(pm@source@pp)/(eta-1j*gap)[:,None,None])

def current(spec,field=0.,direction='xy'):
    if not np.isfinite(field) or direction not in ('xy','yx'):
        raise ValueError('invalid field or current component')
    kx,ky,w=quadrature(spec.nr,spec.nphi,spec.cutoff)
    en,pp,pm,rho,source=_blocks(spec,kx,ky,field)
    if direction=='yx':
        qy=ky+field;rr=np.hypot(kx,qy)
        f=1/(1+np.exp(np.clip((en-spec.mu)/spec.temperature,-500,500)))
        fp=-f*(1-f)/spec.temperature
        df=np.stack((spec.beta+kx/rr,spec.beta-kx/rr),axis=1)
        nder=(kx*qy/rr**3)[:,None,None]*SX+(qy*qy/rr**3)[:,None,None]*SY
        source=(fp[:,0]*df[:,0])[:,None,None]*pp+(fp[:,1]*df[:,1])[:,None,None]*pm+(f[:,0]-f[:,1])[:,None,None]*nder*.5
        v=-SX
    else:v=spec.beta*I+SY
    resp=_inverse(spec,en,pp,pm,source)
    return float(np.real(np.sum(w*np.trace(v@resp,axis1=-2,axis2=-1))))

def decompose(spec,h=1e-3,direction='xy'):
    """d rho_E(B)/dB = L0^-1(d_B D_y f(H_B)) + L0^-1(-i[M,rho_E])."""
    if h<=0 or h>.1 or direction not in ('xy','yx'):raise ValueError('invalid step/component')
    if direction!='xy':raise ValueError('yx split not implemented')
    kx,ky,w=quadrature(spec.nr,spec.nphi,spec.cutoff)
    en,pp,pm,eq,source=_blocks(spec,kx,ky,0.)
    sp=_blocks(spec,kx,ky,h)[4];sm=_blocks(spec,kx,ky,-h)[4]
    rs=_inverse(spec,en,pp,pm,(sp-sm)/(2*h))
    rE=_inverse(spec,en,pp,pm,source)
    M=-SX
    rc=_inverse(spec,en,pp,pm,-1j*(M@rE-rE@M))
    v=spec.beta*I+SY
    integrate=lambda a:float(np.real(np.sum(w*np.trace(v@a,axis1=-2,axis2=-1))))
    direct=(current(spec,h)-current(spec,-h))/(2*h)
    source_val=integrate(rs);comm_val=integrate(rc)
    return dict(protocol='pre_equilibrated_scalar_rta_stationary_E_linear_B_derivative',
      scientific_status='PARTIAL',source_prepare=source_val,magnetic_commutator=comm_val,
      complete_split=source_val+comm_val,full_finite_difference=direct,
      decomposition_residual=source_val+comm_val-direct,eta=spec.eta,beta=spec.beta,
      h=h,nr=spec.nr,nphi=spec.nphi,cutoff=spec.cutoff,
      warnings=['Full response is shift-invariant for fixed-tilt infinite continuum',
                'This declared collision target is not necessarily the source-paper RTA'])

def _b_response_unprepared(spec,kx,ky):
    en,pp,pm,rho,source=_blocks(spec,kx,ky,0.)
    M=-SX
    return _inverse(spec,en,pp,pm,-1j*(M@rho-rho@M))

def unprepared_two_field(spec,step=2e-4):
    """Frozen H0 RTA target: rho_B=L0^-1(-i[M,rho0]); both EB orderings."""
    if not 0<step<.01:raise ValueError('positive small step required')
    kx,ky,w=quadrature(spec.nr,spec.nphi,spec.cutoff)
    en,pp,pm,eq,source=_blocks(spec,kx,ky,0.)
    rE=_inverse(spec,en,pp,pm,source)
    rp=_b_response_unprepared(spec,kx,ky+step)
    rm=_b_response_unprepared(spec,kx,ky-step)
    M=-SX
    be=_inverse(spec,en,pp,pm,(rp-rm)/(2*step))
    eb=_inverse(spec,en,pp,pm,-1j*(M@rE-rE@M))
    v=spec.beta*I+SY
    integrate=lambda a:float(np.real(np.sum(w*np.trace(v@a,axis1=-2,axis2=-1))))
    return dict(protocol='frozen_H0_RTA_two_orderings_field_unprepared',
       scientific_status='PARTIAL',b_then_e=integrate(be),e_then_b=integrate(eb),
       total=integrate(be+eb),eta=spec.eta,step=step,
       limitation='Collision target f(H0) frozen; NOT identified with author SM prescription')

def audit(spec=Spec(),h=.001):
    d=decompose(spec,h)
    xplus=current(spec,h,'xy');xminus=current(spec,-h,'xy')
    yplus=current(spec,h,'yx');yminus=current(spec,-h,'yx')
    d.update(sigma_xy_zeroB=current(spec,0,'xy'),sigma_xy_plusB=xplus,
       sigma_xy_minusB=xminus,sigma_yx_zeroB=current(spec,0,'yx'),
       sigma_yx_plusB=yplus,sigma_yx_minusB=yminus,
       symmetric_cross_derivative=(xplus-xminus+yplus-yminus)/(4*h),
       antisymmetric_cross_derivative=(xplus-xminus-yplus+yminus)/(4*h))
    return d
