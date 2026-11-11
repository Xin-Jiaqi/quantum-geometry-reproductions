"""Shape-aware local selection rule. Does not certify a material or global optimum."""
from __future__ import annotations
import numpy as np

OPTIMA = {
    'gaussian': {'zeta': 2.8412622582895133, 'xi': -.07023878640573547, 'phi': 3.222617743244747},
    'lorentzian': {'zeta': 3.588994310424154, 'xi': -.03176316473924283, 'phi': 2.6648606624504088},
}

def predict(critical: dict, width: float, shape: str='gaussian', *, fwhm: bool=False) -> dict:
    """Predict a control parameter using a quadratic-quartic local energy surface.

    Required critical keys: mass (= generic control lambda_c), A>0, B>0,
    a_prime = da/dlambda != 0, kernel = smooth nonzero K_c.
    width: Gaussian sigma / Lorentzian HWHM unless fwhm=True.
    The returned amplitude excludes regular background and physical SI prefactors.
    """
    if shape not in OPTIMA:raise ValueError('Supported shapes: gaussian, lorentzian')
    if not np.isfinite(width) or width<=0:raise ValueError('Width must be finite and positive')
    values=[critical[k] for k in ('mass','A','B','a_prime','kernel')]
    if not np.isfinite(values).all():raise ValueError('Finite local data required')
    mc,A,B,slope,K=values
    if A<=0 or B<=0:raise ValueError('Requires Delta=Delta_T-A qh^2+a qs^2-B qs^4, A,B>0')
    if abs(slope)<1e-14:raise ValueError('Transverse control crossing is absent; square-root rule does not apply')
    if abs(K)<1e-14:raise ValueError('Nonzero leading response kernel required; JDOS alone is insufficient')
    factor=2*np.sqrt(2*np.log(2)) if shape=='gaussian' else 2.
    scale=width/factor if fwhm else width
    const=OPTIMA[shape];a=const['zeta']*np.sqrt(B*scale)
    return {'shape':shape,'scale':float(scale),'FWHM':float(factor*scale),
            'predicted_control':float(mc+a/slope),'predicted_soft_coefficient':float(a),
            'photon_offset_from_saddle':float(const['xi']*scale),
            'critical_energy_separation':float(const['zeta']**2*scale/4),
            'local_signed_response':float(K*scale**(-.25)/(np.sqrt(A)*B**.25*(2*np.pi)**2)*const['phi']),
            'scope':'Local asymptotic, fixed smooth nonzero kernel, one critical neighborhood; no material/global-optimum certification.'}
