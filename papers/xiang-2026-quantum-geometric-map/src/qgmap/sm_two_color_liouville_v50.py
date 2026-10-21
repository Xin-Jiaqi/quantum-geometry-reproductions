"""v50: direct Xiang et al. original SM two-colour Liouville Eq4-10 & Eq20 pole.

Two INDEPENDENT positive adiabatic regulators eta_a and eta_b and drives
omega_a, omega_b; rho_EB includes both external-field orderings.
Spin perturbation is Bbar*sigma_x, Bbar=-g*muB*physical_B. No uniform
finite scattering/bath is assumed. Current is per single E*Bbar Fourier
amplitude, with physical B sign converted explicitly. All sources are
gauge-independent spin-basis matrices; k derivatives are full finite differences.
Original Eq20 and printed Eq31 signs remain in tension; not silently fixed.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
import numpy as np
from .sm_literal_spin_phe_v49 import SourceConfig,model as sm_literal,SX,SY,ID
from .dirac_mixed_rta_v46 import quadrature

@dataclass(frozen=True)
class Tone:
    omega_a:float=0.
    omega_b:float=0.
    eta_a:float=.001
    eta_b:float=.002
    dk:float=2e-5
    def __post_init__(self):
        if not np.isfinite([self.omega_a,self.omega_b,self.eta_a,self.eta_b,self.dk]).all():
            raise ValueError('finite two-colour parameters')
        if self.eta_a<=0 or self.eta_b<=0 or not 0<self.dk<.01:
            raise ValueError('positive independent adiabatic poles and finite difference required')

def _matrices(x,y,c:SourceConfig,helicity):
    if helicity not in ('main','sm'):raise ValueError('main or sm helicity required')
    sg=1 if helicity=='main' else -1
    H=c.beta*x[:,None,None]*ID+sg*(x[:,None,None]*SY-y[:,None,None]*SX)
    e,u=np.linalg.eigh(H)
    ud=u.conj().swapaxes(-1,-2)
    f=1/(1+np.exp(np.clip((e-c.mu)/c.temperature,-500,500)))
    rho=(u*f[:,None,:])@ud
    return H,e,u,ud,rho

def _resolvent(h,e,u,ud,source,z):
    """Exact inverse of L_z=z+i[H,.] via 2-band eigensystem."""
    if not np.isfinite([np.real(z),np.imag(z)]).all() or np.real(z)<=0:
        raise ValueError('positive regulator required')
    tmp=ud@source@u
    de=e[:,:,None]-e[:,None,:]
    return u@(tmp/(z+1j*de))@ud

def _rho0_source(x,y,c,helicity,dk,direction):
    if direction=='y':xp=x;yp=y+dk;xm=x;ym=y-dk
    elif direction=='x':xp=x+dk;yp=y;xm=x-dk;ym=y
    else:raise ValueError('direction must be x or y')
    return (_matrices(xp,yp,c,helicity)[-1]-_matrices(xm,ym,c,helicity)[-1])/(2*dk)

def _rhoB(x,y,c,helicity,zb):
    H,e,u,ud,rho=_matrices(x,y,c,helicity)
    return _resolvent(H,e,u,ud,-1j*(SX@rho-rho@SX),zb)

def _rhoB_derivative(x,y,c,helicity,tone,direction,zb):
    dk=tone.dk
    if direction=='y':xp=x;yp=y+dk;xm=x;ym=y-dk
    else:xp=x+dk;yp=y;xm=x-dk;ym=y
    return (_rhoB(xp,yp,c,helicity,zb)-_rhoB(xm,ym,c,helicity,zb))/(2*dk)

def two_color(c=SourceConfig(nr=56,nphi=80),tone=Tone(),
              helicity='main',direction='xy',pole_only=False):
    """Direct Eq9 original two-colour Liouville in spin basis.

    rho_E=L(za)^-1 d_k rho0, rho_B=L(zb)^-1(-i[sigma_x,rho0])
    rho_EB=L(za+zb)^-1(d_k rho_B-i[sigma_x,rho_E]).
    Output is per physical B=-Bbar, using source angular/radial quadrature.
    """
    if helicity not in ('main','sm') or direction not in ('xy','yx'):
        raise ValueError('unsupported helicity/direction')
    x,y,w=quadrature(c.nr,c.nphi,c.cutoff)
    h,e,u,ud,rho=_matrices(x,y,c,helicity)
    bdir='y' if direction=='xy' else 'x'
    vx=c.beta*ID+(SY if helicity=='main' else -SY)
    vy=(-SX if helicity=='main' else SX)
    v=vx if direction=='xy' else vy
    za=tone.eta_a-1j*tone.omega_a
    zb=tone.eta_b-1j*tone.omega_b
    zsum=za+zb
    D0=_rho0_source(x,y,c,helicity,tone.dk,bdir)
    if pole_only:
        if tone.omega_a!=0:raise ValueError('electric pole requires omega_a=0')
        fn=1/(1+np.exp(np.clip((e-c.mu)/c.temperature,-500,500)))
        fp=-fn*(1-fn)/c.temperature
        vel=ud@(vy if bdir=='y' else vx)@u
        vdiag=np.real(np.diagonal(vel,axis1=-2,axis2=-1))
        Ddiag=(u*(fp*vdiag)[:,None,:])@ud
        magnetic=-1j*(SX@Ddiag-Ddiag@SX)
        result=_resolvent(h,e,u,ud,magnetic,zsum)
        current=-np.sum(w*np.trace(v@result,axis1=-2,axis2=-1))
        return {'protocol':'SM_Eq20_electric_intraband_pole_at_finite_etas_v50',
                'helicity':helicity,'direction':direction,'tone':asdict(tone),
                'eta_a_times_physical_B_current_residue':complex(current),
                'scientific_status':'PARTIAL'}
    rhoE=_resolvent(h,e,u,ud,D0,za)
    db=_rhoB_derivative(x,y,c,helicity,tone,bdir,zb)
    b_then_e=_resolvent(h,e,u,ud,db,zsum)
    e_then_b=_resolvent(h,e,u,ud,-1j*(SX@rhoE-rhoE@SX),zsum)
    integrate=lambda a:complex(-np.sum(w*np.trace(v@a,axis1=-2,axis2=-1)))
    be=integrate(b_then_e);eb=integrate(e_then_b)
    return {'protocol':'authentic_SM_two_frequency_Liouville_v50',
        'scientific_status':'PARTIAL',
        'source':'SM Eq(4)-(10) both orderings, positive eta_a/eta_b; no all-sector RTA',
        'helicity':helicity,'direction':direction,'tone':asdict(tone),
        'physical_B_b_then_e':be,'physical_B_e_then_b':eb,
        'physical_B_total':be+eb,
        'eta_a_times_b_then_e':tone.eta_a*be,
        'eta_a_times_e_then_b':tone.eta_a*eb,
        'eta_a_times_total':tone.eta_a*(be+eb),
        'limitation':'Unrelaxed two-colour adiabatic protocol. Tau introduced only AFTER selected electric intraband pole; not an equilibrium bath.'}

def run_limit(c=SourceConfig(nr=72,nphi=96),
              eta_as=(.01,.005,.0025,.00125),eta_b=.001,omega_b=0.,
              dk=2e-5,helicity='main',direction='xy'):
    if not eta_as or eta_b<=0:raise ValueError('invalid pole scan')
    rows=[]
    for ea in eta_as:
        t=Tone(eta_a=float(ea),eta_b=float(eta_b),omega_b=omega_b,dk=dk)
        ans=two_color(c,t,helicity,direction)
        pole=two_color(c,t,helicity,direction,pole_only=True)
        rows.append({'eta_a':float(ea),
          'eta_a_times_total_real':ans['eta_a_times_total'].real,
          'eta_a_times_total_imag':ans['eta_a_times_total'].imag,
          'eta_a_times_e_then_b_real':ans['eta_a_times_e_then_b'].real,
          'eta_a_times_b_then_e_real':ans['eta_a_times_b_then_e'].real,
          'intraband_pole_residue_real':pole['eta_a_times_physical_B_current_residue'].real,
          'intraband_pole_residue_imag':pole['eta_a_times_physical_B_current_residue'].imag})
    literal=sm_literal(c,helicity)[direction]
    return {'protocol':'v50_source_SM_Eq9_two_frequency_to_Eq20_pole_and_Eq31',
      'scientific_status':'PARTIAL','config':asdict(c),'eta_b':float(eta_b),
      'omega_b':float(omega_b),'eta_a_scan':rows,
      'helicity':helicity,'direction':direction,
      'SM_literal_Eq20_physical_B':literal['SM_Eq20_current_physical_B'],
      'SM_printed_Eq31_physical_B':literal['SM_printed_Eq31_current_physical_B'],
      'warning':'Eq20/31 original printed signs remain opposite, not patched'}
