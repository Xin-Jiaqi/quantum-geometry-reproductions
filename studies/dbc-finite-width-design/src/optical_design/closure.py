"""Pre-DFT theory closure: normal forms, skew unfolding, error and input gates.

New independent analysis built on the source-independent round-4 local kernel.
The leading theory assumes a smooth optical gap and kernel, not a band degeneracy.
No line width is interpreted as a microscopic collision rate.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from .lineshapes import profile_x


@dataclass(frozen=True)
class NormalForm:
    """Delta = E0 - A*x² + a*y² + C*x*y² - B*y⁴ + b*y.

    A hard-coordinate shear eliminates C and changes the leading quartic
    coefficient. Cubic soft terms must already have been removed at the cusp.
    """
    A: float
    B: float
    a_prime: float
    kernel: float
    control_c: float = 0.
    C: float = 0.

    def __post_init__(self):
        if not np.isfinite(list(vars(self).values())).all():
            raise ValueError('Finite normal-form coefficients required')
        if self.A <= 0 or self.effective_B <= 0:
            raise ValueError('Requires A>0 and B-C^2/(4A)>0; different normal form otherwise')
        if self.a_prime == 0:
            raise ValueError('The control must cross the soft curvature transversely')
        if self.kernel == 0:
            raise ValueError('Nonzero smooth leading kernel required')

    @property
    def effective_B(self):
        return self.B-self.C*self.C/(4*self.A)

    def dimensionless(self, width, control, soft_tilt=0.):
        if not np.isfinite([width,control,soft_tilt]).all() or width <= 0:
            raise ValueError('Positive finite width and finite controls required')
        B=self.effective_B
        return {'zeta': self.a_prime*(control-self.control_c)/np.sqrt(B*width),
                'nu': soft_tilt/(B**.25*width**.75),
                'hard_scale': np.sqrt(width/self.A), 'soft_scale': (width/B)**.25}


def cusp_profile(xi, zeta, nu=0., shape='gaussian'):
    """Integral over both soft-coordinate signs, independent adaptive quadrature.

    nu multiplies Y in -X²+zeta*Y²-Y⁴+nu*Y-xi. nu=0 recovers round 4.
    Supports the declared finite search envelope, not arbitrary huge arguments.
    """
    if shape not in ('gaussian','lorentzian'):
        raise ValueError('Unknown line shape')
    if not np.isfinite([xi,zeta,nu]).all() or max(abs(xi),abs(zeta))>6 or abs(nu)>2:
        raise ValueError('Local envelope: |xi|,|zeta| <= 6 and |nu| <= 2')
    def f(y):
        d=zeta*y*y-y**4-xi
        return float(profile_x(d+nu*y,shape)+profile_x(d-nu*y,shape))
    val,err=quad(f,0,np.inf,epsabs=2e-9,epsrel=2e-9,limit=250)
    return float(val),float(err)


def skew_optimum(nu=0., shape='gaussian'):
    """Candidate optimum in an explicit compact envelope; no global theorem."""
    if not np.isfinite(nu) or abs(nu)>2:
        raise ValueError('Validated skew envelope is |nu| <= 2')
    def f(z):return -cusp_profile(z[0],z[1],nu,shape)[0]
    sols=[minimize(f,p,bounds=[(-4.,6.),(0.,6.)],method='Nelder-Mead',
                  options={'xatol':3e-6,'fatol':2e-9,'maxiter':180})
          for p in [(-.07,2.84),(-1.3,0.),(.6,4.)]]
    good=[s for s in sols if s.success]
    if not good:raise RuntimeError('Skew optimization did not converge')
    s=min(good,key=lambda q:q.fun)
    return dict(nu=float(nu),shape=shape,xi=float(s.x[0]),zeta=float(s.x[1]),
                peak=float(-s.fun),bounds=[[-4.,6.],[0.,6.]],
                multistart_peaks=[float(-q.fun) for q in good])


def local_error_bound(width, area, kernel_error, energy_error, approx_kernel_max,
                      shape='gaussian', omitted_absolute=0.):
    """Conservative absolute integral bound on a finite patch, including measure.

    area means integral d^2k/(2*pi)^2 over that patch. Inputs must be genuine
    uniform sup bounds; sample maxima alone do not certify continuum bounds.
    omitted_absolute bounds all omitted background/tails separately.
    """
    vals=[width,area,kernel_error,energy_error,approx_kernel_max,omitted_absolute]
    if not np.isfinite(vals).all() or width<=0 or min(vals[1:])<0:
        raise ValueError('Finite nonnegative bounds and positive width required')
    if shape=='gaussian':
        l0=1/np.sqrt(2*np.pi); l1=np.exp(-.5)/np.sqrt(2*np.pi)
    elif shape=='lorentzian':
        l0=1/np.pi; l1=9/(8*np.pi*np.sqrt(3))
    else:raise ValueError('Unknown line shape')
    return float(area*(l0*kernel_error/width+l1*approx_kernel_max*energy_error/width**2)
                 +omitted_absolute)


def optimizer_radius(gradient_error, concavity):
    """Strong-concavity conditional radius epsilon_gradient/mu in declared units."""
    if not np.isfinite([gradient_error,concavity]).all() or gradient_error<0 or concavity<=0:
        raise ValueError('Need nonnegative gradient error and positive concavity bound')
    return float(gradient_error/concavity)


def candidate_winner(values, absolute_error):
    """Certify only a finite candidate list when all |F-F_tilde| <= error."""
    v=np.asarray(values,float)
    if v.ndim!=1 or len(v)<2 or not np.isfinite(v).all() or not np.isfinite(absolute_error) or absolute_error<0:
        raise ValueError('Finite candidates and a common nonnegative error bound required')
    order=np.argsort(v);gap=float(v[order[-1]]-v[order[-2]])
    return dict(index=int(order[-1]),gap=gap,certified_among_candidates=bool(gap>2*absolute_error))


def geometry_weight_exponent(p=0,q=0):
    """K proportional to q_h^p q_s^q: exponent, provided its signed moment survives.

    Odd moments may vanish by symmetry. A zero kernel has no guaranteed power.
    """
    if type(p) is not int or type(q) is not int or min(p,q)<0:
        raise ValueError('Nonnegative integer orders required')
    return -.25+p/2+q/4


def audit_material_contract(d):
    """Reject missing scientific inputs; does not calculate a material response."""
    required={'energy_unit':'eV','momentum_basis':'cartesian','response':'normal_shift_current',
              'position_operator':'validated_full','kernel_crosscheck':True,
              'occupation':'insulator','dimension':2,'linewidth_role':'spectral_only',
              'tensor_frame':'fixed_cartesian','bands_resolved':True,'local_symmetry_checked':True,
              'competing_features_checked':True,'full_bz_validation':True}
    errors=[f'{k}: expected {v!r}' for k,v in required.items() if d.get(k)!=v]
    for k in ('input_hash','operator_gauge_id','model_version'):
        if not isinstance(d.get(k),str) or not d[k].strip():errors.append(f'{k}: missing provenance')
    for k in ('insulating_gap_eV','width_eV'):
        if not isinstance(d.get(k),(int,float)) or not np.isfinite(d[k]) or d[k]<=0:
            errors.append(f'{k}: positive finite value required')
    if d.get('line_shape') not in ('gaussian','lorentzian'):errors.append('line_shape: unknown')
    expected={'gaussian':'standard_deviation','lorentzian':'HWHM'}.get(d.get('line_shape'))
    if d.get('width_convention')!=expected or expected is None:errors.append('width_convention: inconsistent')
    return {'accepted':not errors,'errors':errors,
            'scope':'Contract check only; input assertions need external evidence, not self-certification.'}


def quadratic_edge_response(photon, edge_energy, hessian, kernel, width,
                            shape='gaussian', *, multiplicity=1):
    """Signed broadened 2D quadratic minimum/maximum contribution.

    Momentum measure is d²k/(2*pi)². Saddles need an explicit finite patch
    or background prescription and are refused by this cutoff-free formula.
    """
    from scipy.special import ndtr
    H=np.asarray(hessian,float)
    if H.shape!=(2,2) or not np.isfinite(H).all() or not np.allclose(H,H.T,rtol=0,atol=1e-12):
        raise ValueError('Finite symmetric 2x2 Hessian required')
    if not np.isfinite([photon,edge_energy,kernel,width]).all() or width<=0:
        raise ValueError('Finite data and positive width required')
    if type(multiplicity) is not int or multiplicity<1:raise ValueError('Positive valley count required')
    ev=np.linalg.eigvalsh(H)
    if np.prod(ev)<=0:raise ValueError('Nondegenerate minimum/maximum required; saddle is cutoff dependent')
    z=(photon-edge_energy)/width
    if shape=='gaussian':F=float(ndtr(z))
    elif shape=='lorentzian':F=float(.5+np.arctan(z)/np.pi)
    else:raise ValueError('Unknown line shape')
    if ev[0]<0:F=1-F
    return float(multiplicity*kernel*F/(2*np.pi*np.sqrt(np.linalg.det(H))))


def phase_space_exponent(powers):
    """Local definite monomial extremum: measure/delta exponent sum(1/p)-1.

    A saddle can require a cutoff; this helper does not classify saddles.
    """
    if not powers or any(type(p) is not int or p<2 or p%2 for p in powers):
        raise ValueError('Positive even monomial orders required')
    return sum(1/p for p in powers)-1


def massive_dirac_integral(w, u, mass):
    """Round-2 type-I leading kernel, infinite linear patch, frozen optical numerator.

    w is detuning, mass is effective half-gap (not necessarily the orbital m).
    Values at the edge use its one-sided limit; a single point has zero measure.
    """
    if not np.isfinite([w,u,mass]).all() or not 0<u<1 or mass<0:
        raise ValueError('Finite detuning, 0<u<1 and nonnegative mass required')
    if abs(w)<mass*np.sqrt(1-u*u) or w==0:return 0.
    return float(2*np.pi/u*np.sign(w)*(1/np.sqrt(1-u*u)-abs(w)/np.sqrt(w*w+u*u*mass*mass)))
