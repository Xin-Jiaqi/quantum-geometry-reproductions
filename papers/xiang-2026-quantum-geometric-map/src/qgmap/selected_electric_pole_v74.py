"""Source-scoped independent audit of SM Eq.(15)/(20) selected electric pole.

This module independently computes the homogeneous-length-gauge Xiang two-colour
Dirac Liouvillian with exact Fermi-matrix derivatives, and decomposes the
E->B insertion into the electric intraband pole and its complement.
It deliberately DOES NOT claim a dc conductivity, physical bath or Eq31 sign.
Primary-source printed conventions are pinned by the existing v49/v50 audit.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss

I=np.eye(2,dtype=complex)
SX=np.array([[0.,1.],[1.,0.]],complex)
SY=np.array([[0.,-1j],[1j,0.]],complex)

@dataclass(frozen=True)
class Spec:
    beta:float=.125
    mu:float=.3
    temperature:float=.03
    nr:int=52
    nphi:int=72
    cutoff:float=2.5
    helicity:str='main'
    direction:str='xy'
    omega_e:float=0.
    omega_b:float=0.
    eta_e:float=.003
    eta_b:float=.001
    def __post_init__(self):
        real=(self.beta,self.mu,self.temperature,self.cutoff,self.omega_e,
              self.omega_b,self.eta_e,self.eta_b)
        if not np.isfinite(real).all() or not 0<self.beta<1 or min(self.temperature,self.cutoff,self.eta_e,self.eta_b)<=0:
            raise ValueError('finite 0<beta<1, T, cutoff, etaE, etaB positive')
        if not isinstance(self.nr,int) or self.nr<16 or not isinstance(self.nphi,int) or self.nphi<24:
            raise ValueError('quadrature nr>=16 nphi>=24')
        if self.helicity not in ('main','sm') or self.direction not in ('xy','yx'):
            raise ValueError('source helicity and electric/current indices required')

@lru_cache(maxsize=48)
def _grid(nr,nphi,cutoff):
    v,w=leggauss(nr);r=(v+1)*cutoff/2;wr=w*cutoff/2
    phi=2*np.pi*(np.arange(nphi)+.357)/nphi
    x=(r[:,None]*np.cos(phi)).ravel()
    y=(r[:,None]*np.sin(phi)).ravel()
    weights=np.broadcast_to((wr*r)[:,None]/(2*np.pi*nphi),(nr,nphi)).ravel()
    return x,y,weights

def _dag(a):return np.swapaxes(a.conj(),-1,-2)

def _comm(a,b):return a@b-b@a

def _data(s):
    x,y,weight=_grid(s.nr,s.nphi,s.cutoff)
    sg=1 if s.helicity=='main' else -1
    h=s.beta*x[:,None,None]*I+sg*(x[:,None,None]*SY-y[:,None,None]*SX)
    vx=s.beta*I+sg*SY;vy=-sg*SX
    drive=vy if s.direction=='xy' else vx
    current=vx if s.direction=='xy' else vy
    e,u=np.linalg.eigh(h)
    ud=_dag(u)
    f=1/(1+np.exp(np.clip((e-s.mu)/s.temperature,-600,600)))
    fp=-f*(1-f)/s.temperature
    rho0=(u*f[:,None,:])@ud
    d=e[:,:,None]-e[:,None,:]
    gap=f[:,:,None]-f[:,None,:]
    midpoint=1/(1+np.exp(np.clip(((e[:,:,None]+e[:,None,:])/2-s.mu)/s.temperature,-600,600)))
    deriv=-midpoint*(1-midpoint)/s.temperature
    a=np.divide(gap,d,out=np.zeros_like(gap),where=np.abs(d)>1e-11)
    a=np.where(abs(d)>1e-11,a,deriv)
    vb=ud@drive@u
    De=(u@(a*vb)@ud)
    diag_v=np.real(np.diagonal(vb,axis1=-2,axis2=-1))
    Ddiag=(u*(fp*diag_v)[:,None,:])@ud
    Doff=De-Ddiag
    return locals()

def _resolv(d,source,z):
    return d['u']@((_dag(d['u'])@source@d['u'])/(z+1j*d['d']))@d['ud']

def _integrate(d,R):
    return complex(-np.sum(d['weight']*np.trace(d['current']@R,axis1=-2,axis2=-1)))

def _current_from_inputs(s):
    d=_data(s)
    ze=s.eta_e-1j*s.omega_e
    zb=s.eta_b-1j*s.omega_b
    zs=ze+zb
    rho0=d['rho0']
    zspin=SX
    re=_resolv(d,d['De'],ze)
    rb=_resolv(d,-1j*_comm(zspin,rho0),zb)
    # The exact k-derivative of the magnetic resolvent:
    # (zb + i ad_H) (D rb) = -i[Z,D rho0] - i[v_drive,rb]
    # No finite-k difference, no reuse of v50 source's numerical derivative.
    db=_resolv(d,-1j*_comm(zspin,d['De'])-1j*_comm(d['drive'],rb),zb)
    be=_resolv(d,db,zs)
    eb=_resolv(d,-1j*_comm(zspin,re),zs)
    electric_diag_residue=_resolv(d,-1j*_comm(zspin,d['Ddiag']),zs)
    e_diag_pole=electric_diag_residue/ze
    eb_complement=_resolv(d,-1j*_comm(zspin,_resolv(d,d['Doff'],ze)),zs)
    vals=dict(b_then_e=_integrate(d,be),e_then_b=_integrate(d,eb),
              electric_diag_pole=_integrate(d,e_diag_pole),
              electric_diag_pole_residue=_integrate(d,electric_diag_residue),
              eb_non_electric_diag=_integrate(d,eb_complement))
    vals['full']=vals['b_then_e']+vals['e_then_b']
    vals['eb_decomposition_error']=abs(vals['e_then_b']-vals['electric_diag_pole']-vals['eb_non_electric_diag'])
    vals['eta_e_scaled_full']=s.eta_e*vals['full']
    vals['eta_e_scaled_b_then_e']=s.eta_e*vals['b_then_e']
    vals['eta_e_scaled_e_then_b']=s.eta_e*vals['e_then_b']
    vals['eta_e_scaled_intraband_e_pole']=s.eta_e*vals['electric_diag_pole']
    vals['output_regulator']=zs
    return vals

def channels(s=Spec()):
    """Independent analytical (k-Fréchet) source two-color decomposition.

    'eb_non_electric_diag' may have its own output/sum-frequency pole when
    both adiabatic regulators vanish; do NOT globally call it 'regular'.
    """
    if not isinstance(s,Spec):raise TypeError('Spec required')
    z=_current_from_inputs(s)
    z['status']='SOURCE_HOMOGENEOUS_LG_SELECTED_E_POLE_AUDIT_NOT_FULL_DC_PHE'
    return z

def paths(base=Spec()):
    """Two distinct regulator approaches, retaining source electrostatics."""
    curves={}
    for name,scale in [('etaB_equals_etaE',1.),('etaB_tenth_etaE',.1),('etaB_ten_etaE',10.)]:
        rows=[]
        for eta in (.008,.004,.002,.001):
            s=replace(base,omega_e=0.,omega_b=0.,eta_e=eta,eta_b=scale*eta)
            v=channels(s)
            rows.append({'eta_e':eta,'eta_b':scale*eta,
                'scaled_full':v['eta_e_scaled_full'],
                'scaled_B_then_E':v['eta_e_scaled_b_then_e'],
                'scaled_E_then_B':v['eta_e_scaled_e_then_b'],
                'scaled_explicit_e_diag_pole':v['eta_e_scaled_intraband_e_pole']})
        curves[name]=rows
    return curves

def finite_frequency_negative_control(base=Spec()):
    out=[]
    for eta in (.008,.004,.002,.001):
        s=replace(base,omega_e=.15,omega_b=.04,eta_e=eta,eta_b=eta)
        v=channels(s)
        out.append((eta,v['eta_e_scaled_full'],v['eta_e_scaled_intraband_e_pole']))
    return out

def literal_eq20(s=Spec()):
    """Independently contract literal Eq20 in this standalone 2-band model.

    Per physical B, tau=1. This is *not* precomputed using `channels`.
    """
    d=_data(s)
    vband=d['ud']@d['current']@d['u']
    spin=d['ud']@SX@d['u']
    de=d['e'][:,:,None]-d['e'][:,None,:]
    fgrad=d['fp']*np.real(np.diagonal(d['ud']@d['drive']@d['u'],axis1=-2,axis2=-1))
    f_nm=fgrad[:,None,:]-fgrad[:,:,None]
    ident=np.eye(2,dtype=bool)
    rho20=np.where(ident,0.,-spin*f_nm/np.where(ident,1,de))
    return complex(-np.sum(d['weight']*np.einsum('...mn,...nm->...',rho20,vband)))

def printed_eq31(s=Spec()):
    """The printed Eq31 opposite branch, not assumed physically correct."""
    return -literal_eq20(s)
