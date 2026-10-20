"""v20: exact finite-system real-frequency charge-continuity Kubo benchmark.

This is deliberately independent of the SCBA/ladder kernel. A *finite* periodic
2D two-spin-orbital tight-binding model (v17's Bloch function) is Fourier-
transformed into its full real-space Hamiltonian; scalar site disorder is added
exactly rather than averaged. Retarded density-current correlators are evaluated
from ALL exact one-particle eigenstates at fixed temperature.

At finite lattice momentum q and real probe frequency omega, continuity follows
from [H, rho_q] = -Q(q) J_q and the equal-time contact commutator:
 (omega+i*delta)*chi(rho_q, B) - Q(q)*chi(J_q,B) = <[rho_q,B]>.
Q(q)=2*sin(q/2) is an explicitly declared lattice-derivative convention;
J_q=-[H,rho_q]/Q(q). This defines the finite-q longitudinal current, not
a hand-imposed velocity or a transverse PHE current. The response is dimensionless.

For B=rho_-q the contact is zero. For B=J_-q it generally is nonzero.
A q=0 uniform scalar-potential perturbation dynamically gives chi(N,N)=0
while the GRAND-CANONICAL thermal compressibility d<N>/dmu is positive: distinct
preparation/order-of-limits protocols. Neither quantity is Xiang et al.'s PHE.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .warped_disorder_v17 import Lattice,h_and_vel


@dataclass(frozen=True)
class ExactSystem:
    L: int = 5
    beta: float = 0.125
    alpha: float = 0.16
    lam: float = 0.6
    disorder_std: float = 0.13
    disorder_seed: int = 2026
    mu: float = 0.55
    temperature: float = 0.11

    def __post_init__(self):
        if self.L < 3 or self.L > 12 or self.disorder_std < 0 or self.temperature <= 0:
            raise ValueError('finite system requires 3<=L<=12, nonnegative disorder, T>0')


def _lattice_positions(L):
    return np.array([(x,y) for x in range(L) for y in range(L)], dtype=float)


def exact_hamiltonian(p:ExactSystem):
    """Independent real-space H and site coordinates, with exact Bloch Fourier reconstruction."""
    L=p.L; r=_lattice_positions(L); Ns=L*L
    k=(2*np.pi/L)*r
    pp=Lattice(n=12,beta=p.beta,alpha=p.alpha,lam=p.lam)
    hk,_,_=h_and_vel(k,pp,b=0.,protocol='onsite')
    fourier=np.exp(1j*r@k.T)/np.sqrt(Ns)
    u=np.kron(fourier,np.eye(2))
    diag=np.zeros((2*Ns,2*Ns),complex)
    for a in range(Ns):
        diag[2*a:2*a+2,2*a:2*a+2]=hk[a]
    clean=u@diag@u.conj().T
    rng=np.random.default_rng(p.disorder_seed)
    site=rng.normal(0.,p.disorder_std,Ns)
    disorder=np.diag(np.repeat(site,2))
    h=clean+disorder
    if not np.allclose(h,h.conj().T,atol=1e-12):
        raise AssertionError('H must be Hermitian')
    return h,r,site


def charge_operators(h,positions,*,direction=1,q=None):
    """Exact finite-q density and longitudinal current, with discrete Q(q)."""
    if direction not in (0,1):raise ValueError('direction is x/y')
    L=round(np.sqrt(len(positions)))
    if q is None:q=2*np.pi/L
    q=float(q)
    if not np.isfinite(q) or abs(q)>np.pi:
        raise ValueError('finite nonzero momentum within Brillouin zone required')
    qfac=2*np.sin(q/2)
    if abs(qfac)<1e-12:raise ValueError('q=0 current undefined; use uniform_density')
    rq=np.diag(np.repeat(np.exp(-1j*q*positions[:,direction]),2))
    rminus=rq.conj().T
    jq=-(h@rq-rq@h)/qfac
    jminus=jq.conj().T
    return rq,rminus,jq,jminus,qfac


def fermi(e,mu,T):
    x=np.clip((np.asarray(e)-mu)/T,-500,500)
    return 1/(1+np.exp(x))


def exact_retarded_chi(e,vecs,A,B,occupations,omega,delta):
    """Finite-system full Lehmann formula, no SCBA, no vertex truncation."""
    if delta<=0 or not np.isfinite(omega): raise ValueError('positive delta, real finite omega')
    ud=vecs.conj().T
    ab=ud@A@vecs
    bb=ud@B@vecs
    de=e[:,None]-e[None,:]
    ff=occupations[:,None]-occupations[None,:]
    return np.sum(ff*ab*bb.T/(omega+de+1j*delta))


def equilibrium_contact(vecs,occupations,A,B):
    density=(vecs*occupations[None,:])@vecs.conj().T
    return np.trace(density@(A@B-B@A))


def kinematic_ward(p:ExactSystem,*,direction=1,omega=0.27,delta=0.035):
    h,r,site=exact_hamiltonian(p)
    e,u=np.linalg.eigh(h)
    f=fermi(e,p.mu,p.temperature)
    rho,rhobar,j,jbar,qfac=charge_operators(h,r,direction=direction)
    def chi(a,b):return exact_retarded_chi(e,u,a,b,f,omega,delta)
    rr=chi(rho,rhobar);jr=chi(j,rhobar)
    rj=chi(rho,jbar);jj=chi(j,jbar)
    c_rr=equilibrium_contact(u,f,rho,rhobar)
    c_rj=equilibrium_contact(u,f,rho,jbar)
    w_rr=(omega+1j*delta)*rr-qfac*jr-c_rr
    w_rj=(omega+1j*delta)*rj-qfac*jj-c_rj
    algebra_res=float(np.max(np.abs(h@rho-rho@h+qfac*j)))
    # Nonzero scattering eigenvalue widths reflect *static disorder* only;
    # delta is a separate retarded response regularizer, not a collision selfenergy.
    return {
        'protocol':'exact_disordered_finite_2D_Kubo_real_omega_lattice_continuity',
        'L':p.L,'direction':'x' if direction==0 else 'y',
        'q':float(2*np.pi/p.L),'Q':float(qfac),
        'omega':float(omega),'delta':float(delta),
        'rho_rho':rr,'j_rho':jr,'rho_j':rj,'j_j':jj,
        'contact_rho_rho':c_rr,'contact_rho_j':c_rj,
        'ward_density_abs':float(abs(w_rr)),
        'ward_current_abs':float(abs(w_rj)),
        'commutator_abs':algebra_res,
        'bare_missing_current_contact_abs':float(abs((omega+1j*delta)*rj-qfac*jj)),
        'min_occupation':float(f.min()),'max_occupation':float(f.max()),
        'max_site_disorder':float(np.max(abs(site))),
    }


def uniform_density_vs_compressibility(p:ExactSystem,*,omega=0.27,delta=0.035):
    h,_,_=exact_hamiltonian(p)
    e,u=np.linalg.eigh(h)
    f=fermi(e,p.mu,p.temperature)
    N=np.eye(len(e))
    dyn=exact_retarded_chi(e,u,N,N,f,omega,delta)
    chi_thermo=float(np.sum(f*(1-f)/p.temperature))
    # independent finite difference derivative to check the equilibrium factor
    step=1e-5
    fd=float((fermi(e,p.mu+step,p.temperature).sum()-fermi(e,p.mu-step,p.temperature).sum())/(2*step))
    return dict(protocol='dynamic_uniform_scalar_vs_equilibrium_grand_canonical_mu_shift',
                omega=omega,delta=delta,dynamic_uniform_density=dyn,
                equilibrium_compressibility=chi_thermo,
                compressibility_finite_difference=fd,
                contact=equilibrium_contact(u,f,N,N),
                warning='Different protocols, neither is original Xiang EB PHE')
