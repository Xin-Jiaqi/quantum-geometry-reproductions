"""v39: QWZ two-band scalar-disorder SCBA energy derivative and CLEAN-REFERENCE Streda-sea probe.

The computed sea kernel is a deliberately specified SCBA-dressed *bare-current*
Green-function trace. Retarded/advanced energy derivatives have the essential
(I - dSigma/dE) insertion fixed by the same scalar-disorder SCBA equation.
It is NOT the conserving, vertex-complete Kubo-Bastin/Streda conductivity,
let alone the original 2nd order EB response. All units are dimensionless,
conductances in e^2/h with e=hbar=1.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import roots_legendre
from .streda_clean_sea_v38 import qwz_model, sea_per_k, I, SX, SY, SZ, PAULI as PAULI3

PAULI=np.array([I,SX,SY,SZ])

@dataclass(frozen=True)
class SCBAConfig:
    mesh:int=21
    mass:float=-1.
    disorder_variance:float=.12
    eta:float=.12
    mixing:float=.55
    tolerance:float=3e-12
    max_iter:int=800
    def __post_init__(self):
        if self.mesh<9 or self.eta<=0 or self.disorder_variance<0 or not 0<self.mixing<=1 or self.tolerance<=0 or self.max_iter<5:
            raise ValueError('invalid SCBA regularization, grid or iteration parameters')
        if min(abs(self.mass-v) for v in (-2.,0.,2.))<1e-8:
            raise ValueError('QWZ phase transition rejected')


def _invg(matrix):
    return np.linalg.inv(matrix)


def solve_self_energy(H, energy, cfg:SCBAConfig, initial=None):
    """Sigma_R=W mean_k (E+i eta-H-Sigma_R)^(-1)."""
    if H.shape!=(cfg.mesh,cfg.mesh,2,2): raise ValueError('mesh mismatch')
    if not np.isfinite(energy): raise ValueError('nonfinite energy')
    w=cfg.disorder_variance
    sigma=np.zeros((2,2),complex) if initial is None else np.array(initial,complex)
    z=energy+1j*cfg.eta
    for step in range(cfg.max_iter):
        g=_invg(z*I-H-sigma)
        raw=w*np.mean(g,axis=(0,1))
        nxt=(1-cfg.mixing)*sigma+cfg.mixing*raw
        if np.max(abs(nxt-sigma))<cfg.tolerance:
            sigma=nxt
            break
        sigma=nxt
    else:raise RuntimeError('SCBA no convergence')
    gr=_invg(z*I-H-sigma)
    residual=float(np.max(abs(sigma-w*gr.mean(axis=(0,1)))))
    if residual>max(cfg.tolerance*12,1e-9):raise RuntimeError('SCBA residual exceeds tolerance')
    return dict(sigma=sigma,G=gr,iterations=step+1,residual=residual)


def derivative_vertex(G,W):
    """Solve Gamma=I+W mean G Gamma G, i.e. Gamma=I-Sigma'(E).

    An independent 4 x 4 Pauli-basis Bethe-Salpeter solve.  The result is
    energy-density RR vertex, not a current-conserving transverse ladder.
    """
    if W<0:raise ValueError('negative variance')
    matrix=np.empty((4,4),complex)
    for j,b in enumerate(PAULI):
        avg=(G@b@G).mean(axis=(0,1))
        matrix[:,j]=np.einsum('aij,ji->a',PAULI,avg)*(.5*W)
    rhs=np.array([1.,0.,0.,0.],complex)
    coeff=np.linalg.solve(np.eye(4)-matrix,rhs)
    gamma=np.einsum('a,aij->ij',coeff,PAULI)
    resid=float(np.max(abs(gamma-I-W*np.mean(G@gamma@G,axis=(0,1)))))
    return gamma,resid


def at_energy(H,vx,vy,energy,cfg,initial=None):
    sol=solve_self_energy(H,energy,cfg,initial)
    gamma,resid=derivative_vertex(sol['G'],cfg.disorder_variance)
    dG=-sol['G']@gamma@sol['G']
    # sigma'(E)=I-Gamma; choose exact complex-conjugate advanced branch.
    GA=sol['G'].conj().swapaxes(-1,-2)
    dGA=dG.conj().swapaxes(-1,-2)
    trace=lambda m:np.trace(m,axis1=-2,axis2=-1)
    sea=trace(vx@sol['G']@vy@dG-vx@dG@vy@sol['G']-vx@GA@vy@dGA+vx@dGA@vy@GA)
    # A linear SCBA selfenergy is NOT energy constant; naive -G^2 is not exact.
    naive=-sol['G']@sol['G']; naiveA=naive.conj().swapaxes(-1,-2)
    sea_naive=trace(vx@sol['G']@vy@naive-vx@naive@vy@sol['G']-vx@GA@vy@naiveA+vx@naiveA@vy@GA)
    return dict(sigma=sol['sigma'],G=sol['G'],gamma=gamma,derivative=dG,
                vertex_residual=resid,scba_residual=sol['residual'],
                spectral_sea=float(np.real(np.mean(sea))/(4*np.pi)),
                naive_sea=float(np.real(np.mean(sea_naive))/(4*np.pi)),
                gamma_norm=float(np.max(abs(gamma-I))), iterations=sol['iterations'])


def derivative_audit(mesh=21,mass=-1.,W=.12,energy=.7,eta=.12,h=2e-4):
    cfg=SCBAConfig(mesh=mesh,mass=mass,disorder_variance=W,eta=eta)
    H,vx,vy=qwz_model(mesh,mass)
    mid=at_energy(H,vx,vy,energy,cfg)
    plus=solve_self_energy(H,energy+h,cfg,mid['sigma'])
    minus=solve_self_energy(H,energy-h,cfg,mid['sigma'])
    direct=(plus['sigma']-minus['sigma'])/(2*h)
    analytic=I-mid['gamma']
    direct_gr=(plus['G']-minus['G'])/(2*h)
    return dict(protocol='SCBA_energy_derivative_RR_density_vertex',
                mesh=mesh,mass=mass,W=W,energy=energy,eta=eta,
                self_energy_residual=mid['scba_residual'],
                derivative_vertex_residual=mid['vertex_residual'],
                derivative_max_abs=float(np.max(abs(direct-analytic))),
                propagator_derivative_max_abs=float(np.max(abs(direct_gr-mid['derivative']))),
                gamma_minus_identity=mid['gamma_norm'],
                sample_sea_kernel=mid['spectral_sea'],
                sample_sea_naive=mid['naive_sea'],
                warnings=['RR density vertex, not a transverse current vertex',
                          'finite external eta is a resolvent regulator'])


def sea_energy_scan(mesh=21,mass=-1.,mu=0.,eta=.12,W=.12,
                    lower=-7.,quadrature=224):
    """Energy integrated bare-current antisymmetric SCBA sea kernel.

    This is a mathematical kernel, not a conserving disordered Hall conductivity.
    Each E solves energy-dependent SCBA anew. The naive comparison uses SAME
    Sigma(E) but incorrectly holds Sigma fixed on energy differentiation.
    """
    if lower>=mu or quadrature<20:raise ValueError('invalid energy integration')
    cfg=SCBAConfig(mesh=mesh,mass=mass,eta=eta,disorder_variance=W)
    H,vx,vy=qwz_model(mesh,mass)
    if W==0:
        from .streda_clean_sea_v38 import bulk_sea
        # the analytic v38 primitive is an independent W=0 check
        clean=bulk_sea(mesh,mass,mu,eta)
    else:clean=None
    nodes,weights=roots_legendre(quadrature)
    energies=.5*(mu-lower)*nodes+.5*(mu+lower)
    wt=.5*(mu-lower)*weights
    full=naive=0.;max_res=0.;max_ward=0.;max_gamma=0.;max_iter=0
    # continuation initial for nearby energies; still require fixed-point convergence
    sigma=None
    for E,w in zip(energies,wt):
        cur=at_energy(H,vx,vy,float(E),cfg,sigma)
        sigma=cur['sigma']
        full+=w*cur['spectral_sea']
        naive+=w*cur['naive_sea']
        max_res=max(max_res,cur['scba_residual'])
        max_ward=max(max_ward,cur['vertex_residual'])
        max_gamma=max(max_gamma,cur['gamma_norm'])
        max_iter=max(max_iter,cur['iterations'])
    # as in v38: sigma(e^2/h)=2pi * BZ_mean integral kernel
    factor=2*np.pi
    return dict(protocol='energy_dependent_SCBA_bare_current_antisymmetric_sea_ONLY',
                mesh=mesh,mass=mass,mu=mu,eta=eta,W=W,quadrature=quadrature,lower=lower,
                dressed_sea=float(factor*full),naive_constant_sigma_derivative=float(factor*naive),
                derivative_correction=float(factor*(full-naive)),
                clean_reference_if_W0=None if clean is None else clean['sea_sigma_over_e2h'],
                max_scba_residual=max_res,max_density_vertex_residual=max_ward,
                maximum_gamma_deviation=max_gamma,max_scba_iterations=max_iter,
                explicit_missing=['transverse current ladder','contact/seagull vertices',
                                  'complete Kubo-Streda/Bastin I term','crossed diagrams',
                                  'real physical dc order of limits','original EB nonlinear response'])


def exact_finite_disorder_sea(L=4,mass=-1.,W=.12,mu=0.,eta=.12,seeds=6):
    """Independent exact finite-torus, static scalar site-disorder Fermi-sea control.

    Full 2L^2 eigenstates (NO SCBA). The velocity is the derivative of the
    clean Bloch Hamiltonian before Fourier rotation; onsite disorder has zero
    derivative by explicit convention. Ensemble averaging is finite sample.
    This exact *bare Green trace* is comparable as an approximation diagnostic,
    not a proof that SCBA gives the same fully vertex-correct conductivity.
    """
    if L<3 or L>9 or seeds<1 or W<0 or eta<=0:raise ValueError('invalid finite torus')
    N=L*L
    kv=-np.pi+2*np.pi*np.arange(L)/L
    x,y=np.meshgrid(kv,kv,indexing="ij")
    d=np.stack((np.sin(x),np.sin(y),mass+np.cos(x)+np.cos(y)))
    dx=np.stack((np.cos(x),np.zeros_like(x),-np.sin(x)))
    dy=np.stack((np.zeros_like(y),np.cos(y),-np.sin(y)))
    H=np.einsum("a...,aij->...ij",d,PAULI3)
    vx=np.einsum("a...,aij->...ij",dx,PAULI3)
    vy=np.einsum("a...,aij->...ij",dy,PAULI3)
    pts=np.stack(np.meshgrid(kv,kv,indexing='ij'),axis=-1).reshape((-1,2))
    rs=np.stack(np.meshgrid(np.arange(L),np.arange(L),indexing='ij'),axis=-1).reshape((-1,2))
    U=np.kron(np.exp(1j*rs@pts.T)/np.sqrt(N),np.eye(2))
    def periodic_transform(mat):
        blocks=np.zeros((2*N,2*N),complex)
        for j,m in enumerate(mat.reshape((N,2,2))):
            blocks[2*j:2*j+2,2*j:2*j+2]=m
        return U@blocks@U.conj().T
    bare=periodic_transform(H);Ax=periodic_transform(vx);Ay=periodic_transform(vy)
    if np.max(abs(bare-bare.conj().T))>1e-11:raise RuntimeError('Hamiltonian non Hermitian')
    outputs=[]
    for seed in range(seeds):
        rng=np.random.default_rng(1947+seed)
        potential=np.repeat(rng.normal(scale=np.sqrt(W),size=N),2)
        h=bare+np.diag(potential)
        e,v=np.linalg.eigh(h)
        xx=v.conj().T@Ax@v;yy=v.conj().T@Ay@v
        # sum every unordered interband pair's exact two-pole primitive
        from .streda_clean_sea_v38 import _primitive
        accum=0.
        for a in range(2*N):
            for b in range(a+1,2*N):
                de=e[b]-e[a]
                if de<1e-10:continue
                integral=_primitive(np.array([e[a],e[b]]),de,mu,eta)
                accum+=-np.imag(xx[a,b]*yy[b,a])*np.imag(integral)/np.pi
        outputs.append(float(2*np.pi*accum/N))
    return dict(protocol='exact_finite_torus_static_disorder_bare_sea_Lehmann',
                L=L,Nbands=2*N,mass=mass,W=W,mu=mu,eta=eta,seeds=seeds,
                mean=float(np.mean(outputs)),standard_deviation=float(np.std(outputs)),
                samples=outputs,warning='finite torus and finite disorder ensemble; does not include disordered Hall vertex/completeness')
