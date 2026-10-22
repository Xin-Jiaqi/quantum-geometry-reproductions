"""Cancellation-conditioned BZ uncertainty diagnostics (Round 21).

This is NOT a validated numerical error bound on a Kubo integral. It tests
empirical agreement among INDEPENDENT finite quadratures and computes the
exact algebraic amplification of a hypothetical weighted L1 integrand error.

Let I = integral F, L = integral |F| and retention r = |I|/L.
If a numerical perturbation satisfies integral |delta F| <= eps * L,
then |delta I|/|I| <= eps/r whenever I!=0. Thus a 1%-relative net
target needs eps <= .01*r. This mathematical statement does not show that
any finite k mesh has attained eps, and does not address Hamiltonian/
eigensolver/Wannier/scattering systematic errors.

The ensemble audit does NOT symmetrize the primitive output or declare
CONVERGED from coincident estimates. Aliasing can fool all finite grids.
"""
from __future__ import annotations
import numpy as np
from .primitive_bz import primitive_bz_response
from .hexagonal_bz import hexagonal_bz_response
from .group_microscopic import project_frequency_tensor

def _channel(x,component):
    x=np.asarray(x,complex)
    if x.shape!=(2,3):raise ValueError("response tensor must have shape (2,3)")
    if len(component)!=2 or component[0] not in (0,1) or component[1] not in (0,1,2):
        raise ValueError("component must be Cartesian (a=0/1,b=0/1/2)")
    return complex(x[component[0],component[1]])

def cancellation_condition(retention, *,desired_relative_net_error=.01):
    """Necessary L1 local-accuracy budget under the stated perturbation model.

    When retention is zero, relative net error is undefined; the required
    L1 relative accuracy tends to zero. This is a necessary condition, NOT
    an achieved accuracy or an estimate of the true quadrature error.
    """
    r=float(retention);target=float(desired_relative_net_error)
    if not np.isfinite([r,target]).all() or not (0 <= r <= 1) or target<=0:
        raise ValueError("retention must lie in [0,1] and target must be >0")
    return dict(retention=r,condition_amplification=1/r if r>0 else None,
                desired_relative_net_error=target,
                required_relative_l1_integrand_accuracy=target*r,
                certified=False,meaning="conditional mathematical inequality only")

def _spread(vals):
    vals=np.asarray(vals,complex)
    if len(vals)<2:raise ValueError("at least two estimates required")
    return float(max(abs(vals[i]-vals[j]) for i in range(len(vals)) for j in range(i)))

def summarize_independent_estimates(
    *,coarse_midpoint,fine_midpoint,hex_reference,primitive_gauss,
    amplitude_floor=1e-12
):
    """Unprojected complex-component disagreement, with no false confidence.

    The 'anchor' is the finite-order hex value, never the unknown exact
    response. The empirical disagreement is NOT a rigorous error upper bound.
    """
    coarse=np.asarray(coarse_midpoint,complex)
    fine=np.asarray(fine_midpoint,complex)
    if coarse.ndim!=1 or fine.ndim!=1 or len(coarse)<2 or len(fine)<2:
        raise ValueError("at least two coarse and fine midpoint shifts required")
    if not np.all(np.isfinite(coarse)) or not np.all(np.isfinite(fine)):
        raise ValueError("nonfinite midpoint response estimates")
    h,g=complex(hex_reference),complex(primitive_gauss)
    if not np.isfinite([h.real,h.imag,g.real,g.imag]).all():
        raise ValueError("nonfinite independent reference")
    if not np.isfinite(amplitude_floor) or amplitude_floor<=0:
        raise ValueError("positive absolute amplitude floor required")
    mean_f=complex(np.mean(fine));mean_c=complex(np.mean(coarse))
    shift_spread=_spread(fine)
    refinement=abs(mean_f-mean_c)
    cross_spread=_spread([h,g,mean_f])
    empirical=max(shift_spread,refinement,cross_spread)
    relative=empirical/max(abs(h),amplitude_floor)
    return dict(
        finite_hex_anchor=h,
        finest_midpoint_shift_mean=mean_f,
        finest_primitive_gauss=g,
        fine_shift_spread=float(shift_spread),
        midpoint_refinement_difference=float(refinement),
        independent_method_spread=float(cross_spread),
        empirical_disagreement=float(empirical),
        empirical_disagreement_relative_to_finite_anchor=float(relative),
        # Design choice: never certify from a finite collection of
        # same-Hamiltonian quadratures, especially near cancellations.
        numerical_convergence_status="UNVERIFIED",
        is_rigorous_error_upper_bound=False
    )

def _unit_allowed(group,parity,component):
    T=np.zeros((2,3),float)
    T[component[0],component[1]]=1.
    P=project_frequency_tensor(group,T,parity)
    return bool(abs(P[component[0],component[1]])>1e-10)

def shifted_midpoint_ensemble(
    group, *, component,parity="even",
    coarse_order=8,fine_order=12,
    shifts=((0.,0.),(.25,-.20),(-.20,.25),(.13,.11)),
    hex_order=5,gauss_order=12,
    z=.071+.004j,mu=.13,temperature=.17
):
    """Independent Cartesian tensor-component gate for a named MPG.

    Primitive estimates are intentionally UNPROJECTED. For selected groups
    and small orders this is a diagnostic, not a converged physical value.
    """
    if parity not in ("even","odd"):raise ValueError("parity must be even or odd")
    if not isinstance(coarse_order,int) or not isinstance(fine_order,int) or not 2<=coarse_order<fine_order:
        raise ValueError("require 2<=coarse_order<fine_order")
    if not isinstance(hex_order,int) or not isinstance(gauss_order,int) or min(hex_order,gauss_order)<2:
        raise ValueError("hex and Gauss orders must be >=2")
    if len(shifts)<2:raise ValueError("use multiple independent shifts")
    for sh in shifts:
        arr=np.asarray(sh,float)
        if arr.shape!=(2,) or not np.all(np.isfinite(arr)) or np.any(abs(arr)>=.5):
            raise ValueError("midpoint shifts must have two finite entries with |shift|<.5")
    component=tuple(int(i) for i in component)
    # Explicitly check the group projection *only* to report whether the
    # chosen component is symmetry forbidden, never modify an output.
    allowed=_unit_allowed(group,parity,component)
    common=dict(z=z,mu=mu,temperature=temperature)
    coarse=[];fine=[]
    for shift in shifts:
        coarse.append(_channel(primitive_bz_response(
            group,order=coarse_order,method="midpoint",shift=shift,**common)[parity],component))
        fine.append(_channel(primitive_bz_response(
            group,order=fine_order,method="midpoint",shift=shift,**common)[parity],component))
    full_hex=hexagonal_bz_response(group,order=hex_order,**common)[parity]
    full_gauss=primitive_bz_response(group,order=gauss_order,method="gauss",**common)[parity]
    h=_channel(full_hex,component);g=_channel(full_gauss,component)
    comparison=summarize_independent_estimates(coarse_midpoint=coarse,
        fine_midpoint=fine,hex_reference=h,primitive_gauss=g)
    return dict(
        group=group,component=component,parity=parity,
        symmetry_allowed=allowed,
        orders=dict(midpoint_coarse=coarse_order,midpoint_fine=fine_order,
                    hexagon=hex_order,primitive_gauss=gauss_order),
        shifts=[tuple(s) for s in shifts],
        coarse_midpoint=coarse,fine_midpoint=fine,
        comparison=comparison,
        # Also record unprojected off-manifold error as independent diagnostics.
        hex_projector_residual=float(np.max(abs(full_hex-
                                     project_frequency_tensor(group,full_hex,parity)))),
        gauss_projector_residual=float(np.max(abs(full_gauss-
                                     project_frequency_tensor(group,full_gauss,parity))))
    )
