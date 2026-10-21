"""v43: matched noncrossing SCBA *surface-I(a)* versus exact finite-disorder samples.

Strictly limited to the antisymmetric REAL I(a) current-current RA surface
term of the QWZ linear charge Hall response; NOT full Streda I+II, NOT a
conserving physical DC Hall coefficient, NOT original PRB nonlinear E*B PHE.

SCBA uses the SAME onsite scalar Gaussian variance W as v41 exact torus,
and the SAME bare retarded eta, mu and v=dh/dA as v42. A single RA ladder
is inserted in one current vertex; two dressed legs would double-count.
"""
from __future__ import annotations
import numpy as np
from .streda_clean_sea_v38 import qwz_model, I
from .scba_streda_v39 import SCBAConfig, solve_self_energy
from .current_vertex_v40 import pauli_vertex
from .exact_hall_v41 import torus
from .exact_streda_v42 import direct_surface_trace


def scba_surface(mesh=35, mass=-1., W=.12, mu=1.1, eta=.12,
                 offset=(0.,0.)):
    """Same W self-energy and *single-side* RA Bethe-Salpeter surface I(a).

    Hall=(xy-yx)/2, dimensionless e²/h. Returns diagnostic only; sea missing.
    """
    if mesh<9 or not isinstance(mesh,int) or W<0 or eta<=0 or not np.isfinite([mass,W,mu,eta,*offset]).all():
        raise ValueError('invalid SCBA regulator, model or mesh')
    H,vx,vy=qwz_model(mesh,mass,offset=offset)
    cfg=SCBAConfig(mesh=mesh,mass=mass,disorder_variance=W,eta=eta)
    sol=solve_self_energy(H,mu,cfg)
    GR=sol['G']; GA=GR.conj().swapaxes(-1,-2)
    cx=pauli_vertex(GR,GA,vx,W)
    cy=pauli_vertex(GR,GA,vy,W)
    bubble=lambda a,gl,b,gr:complex(np.mean(np.trace(a@gl@b@gr,axis1=-2,axis2=-1)))
    bare_xy=bubble(vx,GR,vy,GA);bare_yx=bubble(vy,GR,vx,GA)
    dressed_xy=bubble(vx,GR,cy['dressed'],GA)
    dressed_yx=bubble(vy,GR,cx['dressed'],GA)
    # Independently dress left leg via adjoint RA reciprocity with AR ladder.
    leftx=pauli_vertex(GA,GR,vx,W)
    lefty=pauli_vertex(GA,GR,vy,W)
    left_xy=bubble(leftx['dressed'],GR,vy,GA)
    left_yx=bubble(lefty['dressed'],GR,vx,GA)
    density=pauli_vertex(GR,GA,np.broadcast_to(I,GR.shape),W)
    sigma=sol['sigma']
    density_from_sigma=I-(sigma-sigma.conj().T)/(2j*eta)
    ward=float(np.max(abs(density['dressed']-density_from_sigma)))
    bare=float(np.real((bare_xy-bare_yx)/2))
    one=float(np.real((dressed_xy-dressed_yx)/2))
    other=float(np.real((left_xy-left_yx)/2))
    if max(cx['residual'],cy['residual'],leftx['residual'],lefty['residual'])>1e-7:
        raise RuntimeError('ladder failed')
    return {'protocol':'QWZ_non_crossing_matched_SCBA_linear_Hall_Ia_only',
        'scientific_status':'PARTIAL','mesh':mesh,'mass':mass,'W':W,'mu':mu,'eta':eta,
        'bare_surface_Ia':bare,'one_ladder_surface_Ia':one,
        'one_ladder_left_surface_Ia':other,'single_vertex_effect':one-bare,
        'right_left_ladder_difference':one-other,
        'scba_residual':float(sol['residual']),
        'ra_density_resolvent_ward_residual':ward,
        'max_current_bse_residual':max(cx['residual'],cy['residual']),
        'min_ra_vertex_condition':min(cx['condition'],cy['condition']),
        'max_ra_vertex_condition':max(cx['condition'],cy['condition']),
        'selfenergy_coefficients':{'real_trace':float(np.real(np.trace(sigma))/2),
                                   'imag_trace':float(np.imag(np.trace(sigma))/2)},
        'units':'e^2/h for I(a) only',
        'missing':['Streda II sea with matched conserving vertices',
                   'full conserving disorder average / crossed diagrams',
                   'thermodynamic metallic DC limit and EB nonlinear response']}


def exact_ensemble_surface(L=4,mass=-1.,W=.12,mu=1.1,eta=.12,
                           seeds=tuple(range(1947,1971)),ntheta=7):
    """Independent direct finite-matrix surface trace with twist average.

    The ntheta^2 twist points are NOT independent disorder realizations;
    uncertainty is estimated ONLY across separate scalar-disorder seeds.
    W=0 must use one clean result and reports zero sampling SEM.
    """
    if not isinstance(L,int) or L<3 or L>9 or ntheta<2 or W<0 or eta<=0:
        raise ValueError('invalid exact finite-torus protocol')
    if not np.isfinite([mass,W,mu,eta]).all():raise ValueError('nonfinite parameters')
    if not seeds or len(set(seeds))!=len(seeds):raise ValueError('unique seeds required')
    means=[]
    use_seeds=(seeds[0],) if W==0 else seeds
    for seed in use_seeds:
        accum=0.
        for i in range(ntheta):
            for j in range(ntheta):
                theta=(2*np.pi*i/ntheta,2*np.pi*j/ntheta)
                H,Jx,Jy=torus(L,mass,W,int(seed),theta)
                accum+=direct_surface_trace(H,Jx,Jy,mu,eta)['surface_Ia_hall_odd']
        means.append(float(accum/(ntheta*ntheta)))
    mean=float(np.mean(means))
    std=float(np.std(means,ddof=1)) if len(means)>1 else 0.
    return {'protocol':'finite_disorder_exact_RA_Hall_surface_fixed_twist_ensemble',
        'scientific_status':'PARTIAL','L':L,'mass':mass,'W':W,'mu':mu,'eta':eta,
        'ntheta':ntheta,'number_independent_disorder_seeds':len(means),
        'twist_grid_per_seed':ntheta**2,'exact_surface_Ia_mean':mean,
        'seed_standard_deviation':std,
        'seed_standard_error':std/np.sqrt(len(means)),
        'per_seed_twist_averages':means,
        'limitation':'Finite torus and finite disorder ensemble; twist points not independent samples; not physical Hall total'}


def comparison(L=4,mesh=35,mass=-1.,W=.12,mu=1.1,eta=.12,
               seeds=tuple(range(1947,1971)),ntheta=7):
    s=scba_surface(mesh,mass,W,mu,eta)
    ex=exact_ensemble_surface(L,mass,W,mu,eta,seeds,ntheta)
    target=ex['exact_surface_Ia_mean']
    dbare=abs(s['bare_surface_Ia']-target)
    dressed=abs(s['one_ladder_surface_Ia']-target)
    return {'protocol':'same_variance_SCBA_RA_ladder_vs_exact_disorder_Ia_v43',
        'scientific_status':'PARTIAL','scba':s,'exact_finite_ensemble':ex,
        'absolute_error_bare':dbare,'absolute_error_single_ladder':dressed,
        'ladder_improves_this_sample':bool(dressed<dbare),
        'discrepancy_over_seed_sem':(s['one_ladder_surface_Ia']-target)/ex['seed_standard_error']
             if ex['seed_standard_error']>0 else None,
        'warning':'One-loop noncrossing I(a) can deviate from finite-sample exact average; neither alone completes I+II or Xiang mixed EB response'}
