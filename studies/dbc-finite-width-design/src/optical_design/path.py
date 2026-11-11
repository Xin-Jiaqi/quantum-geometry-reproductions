"""Stationarity along an actual, coupled control path at fixed spectral width."""
from __future__ import annotations
import numpy as np


def path_stationarity(phi_value, phi_zeta, phi_nu, *, width, A, B, kernel,
                      a, b, a_prime, b_prime=0., A_prime=0., B_prime=0., kernel_prime=0.):
    """Control derivative after photon-energy optimization (Psi_xi=0).

    Width is fixed along the control path. B means the effective quartic
    coefficient. All supplied slopes must belong to the *same physical path*.
    Away from an energy stationary point the omitted Psi_xi term is required.
    This evaluator does not certify that a root is a maximum.
    """
    values=[phi_value,phi_zeta,phi_nu,width,A,B,kernel,a,b,a_prime,b_prime,A_prime,B_prime,kernel_prime]
    if not np.isfinite(values).all() or min(width,A,B)<=0 or kernel==0:
        raise ValueError('Finite path data, positive scales and nonzero kernel required')
    zeta=a/np.sqrt(B*width);nu=b/(B**.25*width**.75)
    zp=a_prime/np.sqrt(B*width)-zeta*B_prime/(2*B)
    np_=b_prime/(B**.25*width**.75)-nu*B_prime/(4*B)
    logp=kernel_prime/kernel-A_prime/(2*A)-B_prime/(4*B)
    pref=kernel*width**(-.25)/((2*np.pi)**2*np.sqrt(A)*B**.25)
    terms=dict(prefactor_change=logp*phi_value,curvature_change=zp*phi_zeta,skew_change=np_*phi_nu)
    return dict(control_derivative=float(pref*sum(terms.values())),terms=terms,
                zeta=float(zeta),nu=float(nu),zeta_prime=float(zp),nu_prime=float(np_),
                scope='At optimized photon energy, fixed spectral width, stated control path only.')
