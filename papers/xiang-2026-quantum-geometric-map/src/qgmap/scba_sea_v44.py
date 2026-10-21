"""v44: explicit *candidate* SCBA sea II plus v43 surface I(a), QWZ.

This is a quantified noncrossing disorder approximation under a declared Green
regulator, **not** a proven conserving microscopic transport response, nor the
original second-order EB magnetotransport of Xiang et al. PRB (2026).

One-body scalar SCBA Sigma(E) and RR energy Ward Gamma0(E)=I-dSigmaR/dE
are solved at each real-energy quadrature node. Spatial velocity is bare in
sea II, motivated by translational Ward at q=0; this is not a proof that the
full disorder-averaged magnetoelectric Hall protocol is conserving.  We never
silently identify the sum with measured Hall conductivity.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import roots_legendre
from .streda_clean_sea_v38 import qwz_model, I
from .scba_streda_v39 import SCBAConfig, solve_self_energy
from .current_vertex_v40 import pauli_vertex
from .scba_surface_v43 import scba_surface
from .exact_hall_v41 import torus
from .exact_streda_v42 import spectral_streda


def _setup(mesh=21,mass=-1.,W=.12,eta=.12,offset=(0.,0.)):
    if not isinstance(mesh,int) or mesh<9 or not np.isfinite([mass,W,eta,*offset]).all() or W<0 or eta<=0:
        raise ValueError('invalid SCBA sea configuration')
    H,vx,vy=qwz_model(mesh,mass,offset=offset)
    cfg=SCBAConfig(mesh=mesh,mass=mass,disorder_variance=W,eta=eta)
    return H,vx,vy,cfg


def sea_kernel(E,H,vx,vy,cfg,*,finite_difference_h=None):
    """At fixed E: exact RR *energy derivative* through same-W ladder.

    dG/dE=-G Gamma0 G; gamma0=I + W <G Gamma0 G>.
    Independent FD branch resolves Sigma(E+-h) and G(E+-h) from scratch.
    Returns d(sigma_II / (e²/h))/dE (BZ mean), identical convention to v42.
    """
    sol=solve_self_energy(H,float(E),cfg)
    GR=sol['G']; GA=GR.conj().swapaxes(-1,-2)
    den=pauli_vertex(GR,GR,np.broadcast_to(I,GR.shape),cfg.disorder_variance)
    gamma=den['dressed']
    dGR=-GR@gamma@GR; dGA=dGR.conj().swapaxes(-1,-2)
    if finite_difference_h is not None:
        h=float(finite_difference_h)
        if h<=0 or not np.isfinite(h):raise ValueError('positive h required')
        plus=solve_self_energy(H,float(E)+h,cfg,initial=sol['sigma'])
        minus=solve_self_energy(H,float(E)-h,cfg,initial=sol['sigma'])
        approximate=(plus['G']-minus['G'])/(2*h)
        deriv_error=float(np.max(abs(approximate-dGR)))
        sigma_error=float(np.max(abs((plus['sigma']-minus['sigma'])/(2*h)-(I-gamma))))
    else:
        deriv_error=sigma_error=None
    tr=lambda x:np.trace(x,axis1=-2,axis2=-1)
    res=tr(vx@GR@vy@dGR-vx@dGR@vy@GR-vx@GA@vy@dGA+vx@dGA@vy@GA)
    # e²/h normalisation exactly as v42: trace/(2N) for N=L² sites,
    # hence 1/2 BZ mean of 2x2 trace.
    val=float(.5*np.real(np.mean(res)))
    return {'sea_density':val,'gamma0_residual':den['residual'],
            'scba_residual':sol['residual'], 'gamma0_minus_I':float(np.max(abs(gamma-I))),
            'dG_finite_difference_error':deriv_error,
            'sigma_finite_difference_error':sigma_error}


def sea_scan(mesh=21,mass=-1.,W=.12,mu=1.1,eta=.12,offset=(0.,0.),
             lower=-12.,quadrature=420):
    """Gauss-Legendre real-energy integration, SCBA selfconsistent at EACH node.

    Parameters set both the selfenergy regulator eta and the Green functions.
    Convergence with lower energy bound and quadrature must be tested.
    """
    H,vx,vy,cfg=_setup(mesh,mass,W,eta,offset)
    if not np.isfinite([mu,lower]).all() or lower>=mu or quadrature<20:
        raise ValueError('invalid energy range or quadrature')
    nodes,weights=roots_legendre(quadrature)
    es=.5*(mu-lower)*nodes+.5*(mu+lower)
    wt=.5*(mu-lower)*weights
    total=0.;max_density=0.;max_scba=0.;max_gamma=0.
    for E,w in zip(es,wt):
        out=sea_kernel(E,H,vx,vy,cfg)
        total+=w*out['sea_density']
        max_density=max(max_density,out['gamma0_residual'])
        max_scba=max(max_scba,out['scba_residual'])
        max_gamma=max(max_gamma,out['gamma0_minus_I'])
    return {'protocol':'QWZ_energy_dependent_SCBA_sea_II_bare_spatial_current_v44',
            'scientific_status':'PARTIAL', 'mesh':mesh,'mass':mass,'W':W,'mu':mu,
            'eta':eta,'lower':lower,'quadrature':quadrature,'offset':list(offset),
            'sea_II_candidate':float(total), 'max_gamma0_residual':max_density,
            'max_selfenergy_residual':max_scba,'max_gamma0_minus_I':max_gamma,
            'units':'e^2/h',
            'warning':'a specified SCBA sea kernel, NOT proven conserving full disordered Hall'}


def candidate_total(mesh=21,mass=-1.,W=.12,mu=1.1,eta=.12,offset=(0.,0.),
                    lower=-12.,quadrature=420):
    """Candidate SCBA I(a)+II: diagnostic addition with stated missing diagrams."""
    surface=scba_surface(mesh,mass,W,mu,eta,offset)
    sea=sea_scan(mesh,mass,W,mu,eta,offset,lower,quadrature)
    return {'protocol':'QWZ_provisional_non_crossing_SCBA_Ia_plus_II_v44',
            'scientific_status':'PARTIAL','surface_Ia':surface['one_ladder_surface_Ia'],
            'sea_II_candidate':sea['sea_II_candidate'],
            'candidate_total':surface['one_ladder_surface_Ia']+sea['sea_II_candidate'],
            'bare_surface_Ia':surface['bare_surface_Ia'],
            'sea':sea,
            'not_proven':['full conserving disorder averaged Streda','transverse EM Ward/contact',
                          'crossed disorder diagrams','thermodynamic metallic DC limits',
                          'Xiang original nonlinear EB PHE']}


def exact_ensemble_total(L=3,mass=-1.,W=.12,mu=1.1,eta=.12,
                         seeds=(1947,1948,1949,1950),ntheta=3):
    """Fixed finite-seed v42 *full* I(a)+II reference.

    Twist averaged PER seed; sampling uncertainty is calculated across SEEDS.
    No interpretation as thermodynamic disorder average.
    """
    if not isinstance(L,int) or L<3 or L>9 or not isinstance(ntheta,int) or ntheta<2 or W<0 or eta<=0:
        raise ValueError('invalid finite reference')
    if not seeds or len(set(seeds))!=len(seeds):raise ValueError('distinct seeds required')
    if not np.isfinite([mass,W,mu,eta]).all():raise ValueError('nonfinite')
    if W==0:
        raise ValueError('W=0 finite torus may have exact orbital degeneracies; use clean_bloch_pair instead')
    sample=[]
    chosen=seeds
    for seed in chosen:
        acc=np.zeros(3)
        for i in range(ntheta):
            for j in range(ntheta):
                theta=(2*np.pi*i/ntheta,2*np.pi*j/ntheta)
                H,vx,vy=torus(L,mass,W,int(seed),theta)
                z=spectral_streda(H,vx,vy,mu,eta)
                acc+=np.array([z['surface_Ia'],z['sea_II'],z['total']])
        sample.append(acc/(ntheta**2))
    a=np.asarray(sample)
    sem=np.std(a,axis=0,ddof=1)/np.sqrt(len(sample)) if len(sample)>1 else np.zeros(3)
    return {'protocol':'v42_exact_finite_disorder_full_Streda_Ia_II_twist_average',
            'scientific_status':'PARTIAL','L':L,'mass':mass,'W':W,'mu':mu,'eta':eta,
            'ntheta':ntheta,'number_disorder_samples':len(sample),
            'surface_Ia':float(a[:,0].mean()),'sea_II':float(a[:,1].mean()),
            'full_total':float(a[:,2].mean()), 'seed_sem_total':float(sem[2]),
            'seed_sem_surface':float(sem[0]),'seed_sem_sea':float(sem[1]),
            'sample_totals':a[:,2].tolist(),
            'warning':'finite sample average; regulator eta is externally supplied; no infinite-size claim'}


def compare_candidate(mesh=21,L=3,mass=-1.,W=.12,mu=1.1,eta=.12,
                      seeds=(1947,1948,1949,1950),ntheta=3,
                      lower=-12.,quadrature=420):
    approx=candidate_total(mesh,mass,W,mu,eta,lower=lower,quadrature=quadrature)
    exact=(clean_bloch_pair(mesh,mass,mu,eta) if W==0 else
           exact_ensemble_total(L,mass,W,mu,eta,seeds,ntheta))
    return {'protocol':'v44_provisional_SCBA_I_plus_II_vs_v42_exact_finite_ensemble',
            'scientific_status':'PARTIAL','scba_candidate':approx,'finite_exact_reference':exact,
            'surface_discrepancy':approx['surface_Ia']-exact['surface_Ia'],
            'sea_discrepancy':approx['sea_II_candidate']-exact['sea_II'],
            'total_discrepancy':approx['candidate_total']-exact['full_total'],
            'warnings':['candidate total is not proven conserving','exact finite seed not thermodynamic',
                        'separate surface/sea parts are not unique mechanisms']}


def clean_bloch_pair(mesh=21,mass=-1.,mu=1.1,eta=.12,offset=(0.,0.)):
    """Separate analytic CLEAN two-band Streda pair primitive (no SCBA/energy integral).

    This is a periodic Bloch-mesh model check, *not* a disorder-ensemble claim.
    """
    H,vx,vy=qwz_model(mesh,mass,offset=offset)
    e,u=np.linalg.eigh(H)
    ux=u.conj().swapaxes(-1,-2)@vx@u
    uy=u.conj().swapaxes(-1,-2)@vy@u
    amp=np.imag(ux[...,0,1]*uy[...,1,0]); delta=e[...,1]-e[...,0]
    zlo=mu-e[...,0]+1j*eta; zhi=mu-e[...,1]+1j*eta
    primitive=-(1/zlo+1/zhi)/delta+2*(np.log(zlo)-np.log(zhi))/delta**2
    sea=float(-2*np.mean(amp*np.imag(primitive)))
    surf=float(2*np.mean(amp*eta*delta/(abs(zlo)**2*abs(zhi)**2)))
    return {'protocol':'clean_periodic_Bloch_pair_analytic_v44',
            'surface_Ia':surf,'sea_II':sea,'full_total':surf+sea,
            'scientific_status':'PARTIAL','not_disorder_average':True}
