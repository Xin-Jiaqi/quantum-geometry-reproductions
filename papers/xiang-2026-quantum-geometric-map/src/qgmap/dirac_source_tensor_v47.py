"""v47: correct tensor symmetry and independently test source/finite-B response.

The original tilted-Dirac spin PHE is SYMMETRIC in xy and yx. A previous
v46 analytic source used the wrong sign for d_y(kx/r); v47 rechecks the
entire f(H_B) source by finite difference, evaluates both tensor directions
in a matched frozen-H0 RTA, and uses a separately diagonalized finite-B
Liouvillian to validate the mixed-order response.
The authors' original Supplemental Material collision recipe is unavailable,
so a numerical Eq(11) match is a bounded model reproduction, not proof of
the exact author relaxation operator. Scientific status remains PARTIAL.
"""
from __future__ import annotations
import numpy as np
from .dirac_mixed_rta_v46 import (
    Spec,quadrature,_blocks,_inverse,_b_response_unprepared,unprepared_two_field,
    I,SX,SY)
from .dirac_phe_protocol_v45 import printed_eq11_geometry


def derivative_source_check(spec=Spec(nr=16,nphi=32),step=1e-5):
    """Independent matrix finite difference of f(H_B) for the y source."""
    if step<=0 or not np.isfinite(step):raise ValueError('positive step required')
    x,y,_=quadrature(spec.nr,spec.nphi,spec.cutoff)
    analytic=_blocks(spec,x,y,0.)[4]
    fd=(_blocks(spec,x,y+step,0.)[3]-_blocks(spec,x,y-step,0.)[3])/(2*step)
    return {'max_source_y_fd_error':float(np.max(abs(analytic-fd))),
            'max_source_y_fd_value':float(np.max(abs(fd)))}


def unprepared_yx(spec=Spec(),step=3e-5):
    """Magnetic drive on f(H0) target: j_y from E_x, BOTH EB orderings.

    Unlike the analytic y source, the x source is independently differentiated
    from the complete two-band equilibrium matrix. No assumption xy=yx.
    """
    if not np.isfinite(step) or not 0<step<.01:raise ValueError('step')
    x,y,w=quadrature(spec.nr,spec.nphi,spec.cutoff)
    en,pp,pm,rho,_=_blocks(spec,x,y,0.)
    electric=(_blocks(spec,x+step,y,0.)[3]-_blocks(spec,x-step,y,0.)[3])/(2*step)
    rE=_inverse(spec,en,pp,pm,electric)
    bp=_b_response_unprepared(spec,x+step,y)
    bm=_b_response_unprepared(spec,x-step,y)
    be=_inverse(spec,en,pp,pm,(bp-bm)/(2*step))
    M=-SX
    eb=_inverse(spec,en,pp,pm,-1j*(M@rE-rE@M))
    tr=lambda z:float(np.real(np.sum(w*np.trace(-SX@z,axis1=-2,axis2=-1))))
    return {'b_then_e':tr(be),'e_then_b':tr(eb),'total':tr(be+eb),
            'scientific_status':'PARTIAL'}


def tensor(spec=Spec(),step=3e-5):
    """Full xy/yx decomposition. Eq11 comparison ONLY uses symmetric part."""
    xy=unprepared_two_field(spec,step)
    yx=unprepared_yx(spec,step)
    sym=.5*(xy['total']+yx['total'])
    anti=.5*(xy['total']-yx['total'])
    paper=printed_eq11_geometry(spec.beta)*np.tanh(spec.mu/(2*spec.temperature))
    return dict(protocol='corrected_frozen_H0_RTA_symmetric_PHE_v47',
        scientific_status='PARTIAL',beta=spec.beta,eta=spec.eta,mu=spec.mu,
        temperature=spec.temperature,nr=spec.nr,nphi=spec.nphi,step=step,
        xy=xy,yx=yx,symmetric=sym,antisymmetric=anti,
        eta_times_symmetric=spec.eta*sym,eta_times_antisymmetric=spec.eta*anti,
        printed_eq11_geometry=float(paper),
        ratio_to_main_eq11=float(spec.eta*sym/paper) if abs(paper)>1e-13 else None,
        warnings=[
         'v46 electric y source had wrong sign for d_y(kx/r); old frozen data REVOKED',
         'xy alone is NOT the source PHE; compare symmetric (xy+yx)/2',
         'the exact Supplemental Material collision operator is NOT yet checked',
         'Eq8 signed convention and Fig2(b) absolute amplitude remain OPEN'])


def _H(spec,x,y,b):
    return spec.beta*x[:,None,None]*I+x[:,None,None]*SY-(y+b)[:,None,None]*SX


def _density_target(spec,x,y):
    e,U=np.linalg.eigh(_H(spec,x,y,0))
    occ=1/(1+np.exp(np.clip((e-spec.mu)/spec.temperature,-500,500)))
    return (U*occ[:,None,:])@U.conj().swapaxes(-1,-2)


def _direct_superinverse(h,rhs,eta):
    """Second solver: full eigensystem of H_b, not v46 spectral projectors."""
    e,U=np.linalg.eigh(h)
    Ut=U.conj().swapaxes(-1,-2)
    z=(Ut@rhs@U)/(eta+1j*(e[:,:,None]-e[:,None,:]))
    return U@z@Ut


def independent_finite_B_current(spec,b=0.,axis='y',dk=1e-4):
    """Solve frozen-target stationary B density then E directly at finite B.

    rho_B=(eta+i ad_Hb)^-1 eta f(H0);
    rho_E=(eta+i ad_Hb)^-1 d_k rho_B, from a matrix finite difference.
    """
    if axis not in ('x','y') or not np.isfinite([b,dk]).all() or dk<=0:
        raise ValueError('axis/dk')
    x,y,w=quadrature(spec.nr,spec.nphi,spec.cutoff)
    dx=dk if axis=='x' else 0.
    dy=dk if axis=='y' else 0.
    def rb(xx,yy):
        return _direct_superinverse(_H(spec,xx,yy,b),
                    spec.eta*_density_target(spec,xx,yy),spec.eta)
    source=(rb(x+dx,y+dy)-rb(x-dx,y-dy))/(2*dk)
    driven=_direct_superinverse(_H(spec,x,y,b),source,spec.eta)
    v=spec.beta*I+SY if axis=='y' else -SX
    return float(np.real(np.sum(w*np.trace(v@driven,axis1=-2,axis2=-1))))


def independent_mixed_B_derivative(spec=Spec(),axis='y',db=1e-4,dk=1e-4):
    if db<=0:raise ValueError('db')
    return ((independent_finite_B_current(spec,db,axis,dk)-
              independent_finite_B_current(spec,-db,axis,dk))/(2*db))
