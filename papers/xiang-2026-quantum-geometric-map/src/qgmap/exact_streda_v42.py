"""v42 exact finite-torus *linear* Hall Smrcka-Streda decomposition.

Microscopic H, currents and static scalar disorder are imported from v41.
T=0, mu fixed, e=hbar=a=1. Trace/Nsite is a real-space single-sample
trace; eta is a common single-particle Green-function regulator, NOT a
self-consistent scattering rate. This is NOT a conserving SCBA average, and
not original Xiang 2026 mixed E*B second-order magnetotransport.

Smrcka-Streda I(a), I(b), II split can contain a surface-sea overlap: only
the sum is a transport coefficient. Do not equate the individual labels with
separately measurable, unique 'extrinsic'/'intrinsic' physical mechanisms.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad
from .exact_hall_v41 import torus, hall_lehmann


def _input(H,Jx,Jy,mu,eta):
    H=np.asarray(H,dtype=complex)
    Jx=np.asarray(Jx,dtype=complex)
    Jy=np.asarray(Jy,dtype=complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape!=Jx.shape or H.shape!=Jy.shape:
        raise ValueError('Hermitian square H,Jx,Jy of identical shape required')
    sites=H.shape[0]//2
    L=int(round(np.sqrt(sites)))
    if 2*L*L!=H.shape[0] or L<3:
        raise ValueError('two-orbital LxL torus required')
    if not np.isfinite([mu,eta]).all() or eta<=0:
        raise ValueError('finite common eta>0 required')
    if any(not np.isfinite(M).all() or not np.allclose(M,M.conj().T,atol=1e-10)
           for M in (H,Jx,Jy)):
        raise ValueError('finite Hermitian operators required')
    return sites


def _spectral(H,Jx,Jy):
    eig,U=np.linalg.eigh(H)
    x=U.conj().T@Jx@U
    y=U.conj().T@Jy@U
    n,m=np.triu_indices(eig.size,1)
    D=eig[m]-eig[n]
    if np.min(D)<1e-11:
        raise ValueError('nondegenerate full spectrum required for pair primitive')
    P=np.imag(x[n,m]*y[m,n])
    return eig,n,m,D,P


def _sea_primitive(el,eu,mu,eta):
    """Retarded two-pole antiderivative (energy -infinity -> mu)."""
    D=eu-el
    zl=mu-el+1j*eta
    zu=mu-eu+1j*eta
    return -(1/zl+1/zu)/D+2*(np.log(zl)-np.log(zu))/D**2


def spectral_streda(H,Jx,Jy,mu=0.,eta=.04):
    """Pairwise exact Streda I(a) Hall-odd and II; no numerical energy integral.

    sigma/(e²/h), Hall = (xy-yx)/2. I(b) is Hall-even for the
    common Hermitian Hamiltonian and pair of charge current operators.
    """
    ns=_input(H,Jx,Jy,mu,eta)
    eig,n,m,D,A=_spectral(H,Jx,Jy)
    zl=mu-eig[n]+1j*eta
    zu=mu-eig[m]+1j*eta
    prim=_sea_primitive(eig[n],eig[m],mu,eta)
    sea=float(-2/ns*np.sum(A*np.imag(prim)))
    surface=float(2/ns*np.sum(A*eta*D/
        (((mu-eig[n])**2+eta**2)*((mu-eig[m])**2+eta**2))))
    total=surface+sea
    return dict(surface_Ia=surface,sea_II=sea,total=total,
        surface_Ib_hall_odd=0.,eta=eta,mu=mu,
        nearest_level_distance=float(np.min(abs(eig-mu))),
        min_spacing=float(np.min(D)))


def direct_surface_trace(H,Jx,Jy,mu=0.,eta=.04):
    """Independent inverse-matrix I(a), I(b) Hall-odd traces; no eigensolver."""
    ns=_input(H,Jx,Jy,mu,eta)
    GR=np.linalg.inv((mu+1j*eta)*np.eye(H.shape[0])-H)
    GA=GR.conj().T
    tr=np.trace
    Ia_xy=tr(Jx@GR@Jy@GA)/ns
    Ia_yx=tr(Jy@GR@Jx@GA)/ns
    Ib_xy=-tr(Jx@GR@Jy@GR+Jx@GA@Jy@GA)/(2*ns)
    Ib_yx=-tr(Jy@GR@Jx@GR+Jy@GA@Jx@GA)/(2*ns)
    return dict(surface_Ia_hall_odd=float(np.real((Ia_xy-Ia_yx)/2)),
        surface_Ib_hall_odd=float(np.real((Ib_xy-Ib_yx)/2)),
        surface_Ib_hall_odd_complex_abs=float(abs((Ib_xy-Ib_yx)/2)),
        surface_Ia_xy_imag_abs=float(abs(np.imag(Ia_xy))),
        surface_Ia_yx_imag_abs=float(abs(np.imag(Ia_yx))))


def _sea_matrix_integrand(E,H,Jx,Jy,eta):
    """Direct GR/GA and dG/dE; output d sigma_II/(e²/h) / dE."""
    GR=np.linalg.inv((E+1j*eta)*np.eye(H.shape[0])-H)
    GA=GR.conj().T
    dGR=-GR@GR;dGA=-GA@GA
    tr=np.trace
    z=tr(Jx@GR@Jy@dGR-Jx@dGR@Jy@GR-
         Jx@GA@Jy@dGA+Jx@dGA@Jy@GA)
    return float(np.real(z)/(2*(H.shape[0]//2)))


def direct_sea_quadrature(H,Jx,Jy,mu=0.,eta=.1,lower=-35.,epsabs=2e-8):
    """Independent dense Green trace quadrature. NO spectral primitive reused.

    Finite lower cutoff is explicit and must be converged before acceptance.
    """
    _input(H,Jx,Jy,mu,eta)
    if not np.isfinite(lower) or lower>=mu or epsabs<=0:
        raise ValueError('finite lower energy cutoff below mu required')
    poles=np.linalg.eigvalsh(H)
    points=[float(x) for x in poles if lower<x<mu]
    val,err=quad(lambda E:_sea_matrix_integrand(E,H,Jx,Jy,eta),
                 lower,mu,points=points,epsabs=epsabs,epsrel=3e-8,limit=500)
    return dict(sea_II=float(val),integration_error=float(err),lower=lower)


def streda_audit(L=4,mass=-1.,W=.12,seed=1947,theta=(.11,.33),
                 mu=0.,etas=(.16,.05,.01,.003)):
    """One exact torus/sample, full Streda I/II eta scan vs v41 Lehmann DC.

    The last limit is meaningful only if mu is not an eigenlevel.
    """
    H,Jx,Jy=torus(L,mass,W,seed,theta)
    exact=hall_lehmann(H,Jx,Jy,mu)
    rows=[]
    for eta in etas:
        d=spectral_streda(H,Jx,Jy,mu,eta)
        rows.append({'eta':float(eta),'surface_Ia':d['surface_Ia'],
                     'sea_II':d['sea_II'],'total':d['total'],
                     'abs_error_vs_exact_lehmann':abs(d['total']-exact['hall_real'])})
    return dict(protocol='finite_exact_disorder_Smrcka_Streda_linear_Hall_v42',
        scientific_status='PARTIAL',L=L,mass=mass,W=W,seed=seed,
        theta=list(theta),mu=mu,etas=list(etas),
        sigma_units='e^2/h',reference_v41_exact_DC=exact['hall_real'],
        minimum_fermi_gap=exact['minimum_fermi_gap'],
        nearest_level_distance=spectral_streda(H,Jx,Jy,mu,etas[-1])['nearest_level_distance'],
        rows=rows,
        limitations=['finite sample and scalar Green regulator, not SCBA disorder averaging',
          'Streda I/II partition has overlap and need not identify distinct mechanisms',
          'a fixed external eta is not a collision operator or conserving SCBA',
          'the original paper nonlinear EB steady-state DC remains OPEN'])


def ensemble_audit(L=4,mass=-1.,W=.12,seeds=(1947,1948,1949,1950,1951,1952),
                   theta=(.11,.33),mu=1.1,eta=.005):
    """Finite *fixed list* sample mean, no claim of thermodynamic SCBA averaging."""
    if not seeds or len(set(seeds))!=len(seeds):
        raise ValueError('nonempty unique disorder seeds required')
    samples=[]
    for seed in seeds:
        H,Jx,Jy=torus(L,mass,W,int(seed),theta)
        exact=hall_lehmann(H,Jx,Jy,mu)['hall_real']
        st=spectral_streda(H,Jx,Jy,mu,eta)
        samples.append({'seed':int(seed),'exact_lehmann':exact,
                        'surface_Ia':st['surface_Ia'],'sea_II':st['sea_II'],
                        'total':st['total'],
                        'total_minus_lehmann':st['total']-exact,
                        'nearest_level_distance':st['nearest_level_distance']})
    avg=lambda name:float(np.mean([x[name] for x in samples]))
    return dict(protocol='exact_finite_sample_ensemble_Streda_linear_Hall_v42',
                scientific_status='PARTIAL',L=L,mass=mass,W=W,mu=mu,
                theta=list(theta),eta=eta,nsamples=len(seeds),
                samples=samples,mean_exact_lehmann=avg('exact_lehmann'),
                mean_surface_Ia=avg('surface_Ia'),mean_sea_II=avg('sea_II'),
                mean_streda_total=avg('total'),
                std_exact_lehmann=float(np.std([x['exact_lehmann'] for x in samples],ddof=1)) if len(samples)>1 else None,
                mean_total_minus_lehmann=avg('total_minus_lehmann'),
                warnings=['sample average is not infinite-system SCBA or full disorder ensemble',
                          'finite eta and finite level spacing can strongly affect metallic finite-size response',
                          'no separately measurable intrinsic/extrinsic assignment to I/II'])
