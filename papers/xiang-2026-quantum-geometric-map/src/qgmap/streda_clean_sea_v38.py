"""v38: independent CLEAN two-band Kubo-Streda Fermi-sea Hall benchmark.

Important scope: a two-band periodic band insulator, energy-independent scalar
retarded regulator eta>0, no disorder self-energy, no ladder or crossed vertex,
no mixed EB second-order response.  This is the antisymmetric *sea* term, NOT
complete metallic / disordered Kubo-Streda nor paper's DC PHE.

We use the convention Omega_lower=2 Im[(v_x)_{lu}(v_y)_{ul}]/Delta**2,
consistent with the oriented link-plaquette Chern number below. The model units
are e=hbar=1 and sigma_xy/(e^2/h)=C=2*pi*BZ_average(Omega).

For each nondegenerate 2x2 Hamiltonian, the zero-T clean sea trace is
 F(epsilon)=Tr[v_x G^R v_y dG^R - v_x dG^R v_y G^R -(R->A)]/(4pi)
where dG^R/d epsilon=-(G^R)^2 at fixed scalar eta. It is evaluated in two
INDEPENDENT ways: 1) an analytic primitive of the two-pole Green functions,
2) direct real-energy quadrature of the matrix trace. The analytic primitive
contains log(z_l)-log(z_u), with Im z_l,u=eta>0, and zero integration constant
at epsilon -> -infinity.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad

SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.array([[1,0],[0,-1]],complex)
I=np.eye(2,dtype=complex)
PAULI=np.stack([SX,SY,SZ])


def qwz_model(n=61,mass=-1.0,offset=(0.17,0.31)):
    """H=sin(kx) sx+sin(ky) sy+(m+cos(kx)+cos(ky)) sz; analytic velocities.

    Chern with our orientation: -1 for -2<m<0 and +1 for 0<m<2.
    Reject exactly at a bulk gap closing (m=-2,0,+2).
    """
    if n<9: raise ValueError('n>=9 required')
    if np.min(abs(np.asarray([-2.,0.,2.])-mass))<1e-8:
        raise ValueError('QWZ phase transition: gapless at m=-2,0,2')
    a=(-np.pi+2*np.pi*(np.arange(n)+offset[0])/n)
    b=(-np.pi+2*np.pi*(np.arange(n)+offset[1])/n)
    x,y=np.meshgrid(a,b,indexing='ij')
    d=np.stack((np.sin(x),np.sin(y),mass+np.cos(x)+np.cos(y)))
    dx=np.stack((np.cos(x),np.zeros_like(x),-np.sin(x)))
    dy=np.stack((np.zeros_like(x),np.cos(y),-np.sin(y)))
    H=np.einsum('a...,aij->...ij',d,PAULI)
    vx=np.einsum('a...,aij->...ij',dx,PAULI)
    vy=np.einsum('a...,aij->...ij',dy,PAULI)
    return H,vx,vy


def bands(H,vx,vy,minimum_gap=1e-7):
    """Band energies and the complex interband velocity product; no Berry gauge derivatives."""
    H=np.asarray(H,complex)
    vx=np.asarray(vx,complex);vy=np.asarray(vy,complex)
    if H.shape[-2:]!=(2,2) or H.shape!=vx.shape or H.shape!=vy.shape:
        raise ValueError('H/vx/vy must have identical (...,2,2) shapes')
    if not all(np.allclose(z,z.swapaxes(-1,-2).conj(),atol=1e-10) for z in (H,vx,vy)):
        raise ValueError('Hermitian matrices required')
    e,u=np.linalg.eigh(H)
    D=e[...,1]-e[...,0]
    if np.min(D)<minimum_gap:raise ValueError('nondegenerate bands required')
    ud=u.swapaxes(-1,-2).conj()
    wx=ud@vx@u;wy=ud@vy@u
    prod=wx[...,0,1]*wy[...,1,0]
    return e,D,prod,u


def clean_berry(H,vx,vy):
    e,D,prod,_=bands(H,vx,vy)
    lower=2*np.imag(prod)/D**2
    return np.stack((lower,-lower),axis=-1)


def _primitive(e,delta,mu,eta):
    lower=e[...,0];upper=e[...,1]
    z0=mu-lower+1j*eta;z1=mu-upper+1j*eta
    # ∫_{-∞}^mu Δ/(z0² z1²) dε, with Δ=e_upper-e_lower>0
    return -(1/z0+1/z1)/delta + (2/delta**2)*(np.log(z0)-np.log(z1))


def sea_per_k(H,vx,vy,mu=0.,eta=.05):
    """Analytic full sea energy integral in local Berry-curvature units."""
    if eta<=0 or not np.isfinite([mu,eta]).all():raise ValueError('finite eta>0 required')
    e,D,prod,_=bands(H,vx,vy)
    integral=_primitive(e,D,mu,eta)
    # (F^R-F^A)/(4π) integrated: -(Im vx_lu vy_ul)*(Im primitive)/π.
    return -np.imag(prod)*np.imag(integral)/np.pi



def boundary_pole_diagnostic(H,vx,vy,mu=0.,eta=.05):
    """Exact compensating boundary/pole term of our clean sea primitive.

    This is NOT asserted to equal the full Streda Fermi-surface conductivity:
    its role is to isolate the rational (1/z) terms from the sea integral.
    A conserving disordered Fermi-surface sector requires its own Kubo theory.
    """
    if eta<=0:raise ValueError('positive eta')
    e,D,A,_=bands(H,vx,vy)
    lor=(1/((mu-e[...,0])**2+eta**2)+1/((mu-e[...,1])**2+eta**2))
    return np.imag(A)*eta*lor/(np.pi*D)


def log_occupancy_diagnostic(H,vx,vy,mu=0.,eta=.05):
    """Sea plus explicit compensating boundary: broadened Berry occupation.

    This algebraic decomposition is NOT a derivation of a physical Streda I term.
    """
    return sea_per_k(H,vx,vy,mu,eta)+boundary_pole_diagnostic(H,vx,vy,mu,eta)


def metallic_pole_audit(n=101,mass=-1.,mu=1.1,eta=.02):
    H,vx,vy=qwz_model(n,mass)
    sea=float(2*np.pi*np.mean(sea_per_k(H,vx,vy,mu,eta)))
    boundary=float(2*np.pi*np.mean(boundary_pole_diagnostic(H,vx,vy,mu,eta)))
    log_only=float(2*np.pi*np.mean(log_occupancy_diagnostic(H,vx,vy,mu,eta)))
    occ=bulk_sea(n,mass,mu,eta)['occupied_berry_chern']
    return dict(sea_only=sea,compensating_boundary_pole=boundary,
       log_occupancy=log_only,occupied_berry_clean=occ,
       sea_plus_boundary_residual=sea+boundary-log_only,
       warning='Boundary pole is an exact algebraic diagnostic; no complete disordered Streda I vertex')

def bulk_sea(n=61,mass=-1.,mu=0.,eta=.05):
    H,vx,vy=qwz_model(n,mass)
    f=sea_per_k(H,vx,vy,mu,eta)
    bc=clean_berry(H,vx,vy)
    e=np.linalg.eigvalsh(H)
    full_berry=np.sum(bc*(e<mu),axis=-1)
    return dict(model='QWZ_clean_two_band',method='bare_clean_Fermi_sea_only',
        mesh=n,mass=mass,mu=mu,eta=eta,
        sea_sigma_over_e2h=float(2*np.pi*np.mean(f)),
        occupied_berry_chern=float(2*np.pi*np.mean(full_berry)),
        minimum_direct_gap=float(np.min(e[...,1]-e[...,0])),
        abs_difference_to_clean_limit=float(abs(2*np.pi*np.mean(f)-2*np.pi*np.mean(full_berry))),
        warning='finite eta no scattering self-energy/vertex; not complete disordered Kubo-Streda or EB response')


def sea_matrix_integrand(E,H,vx,vy,eta):
    """Direct trace formula, no eigendecomposition; scalar eta, dG=-G²."""
    G=np.linalg.inv((E+1j*eta)*I-H)
    GA=G.conj().swapaxes(-1,-2)
    dG=-G@G;dGA=-GA@GA
    trace=lambda z:np.trace(z,axis1=-2,axis2=-1)
    f=trace(vx@G@vy@dG-vx@dG@vy@G-vx@GA@vy@dGA+vx@dGA@vy@GA)
    return float(np.real(f)/(4*np.pi))


def direct_energy_sea(H,vx,vy,mu=0.,eta=.12,lower_cut=-120.):
    """Slow, independent quadrature used for single k-point acceptance tests."""
    H=np.asarray(H);vx=np.asarray(vx);vy=np.asarray(vy)
    if H.shape!=(2,2):raise ValueError('only a single 2x2 k point')
    if eta<=0:raise ValueError('eta positive')
    if lower_cut>=mu:raise ValueError('integration bound below mu')
    e=np.linalg.eigvalsh(H)
    pts=[float(z) for z in e if lower_cut<z<mu]
    val,err=quad(lambda E: sea_matrix_integrand(E,H,vx,vy,eta),lower_cut,mu,
                 points=pts,epsabs=2e-11,epsrel=3e-11,limit=300)
    return val,err


def fhs_chern(H,band=0,phases=None):
    """Fukui-Hatsugai-Suzuki discrete gauge-invariant link Chern number.

    Works on a periodic 2D mesh, with complex U(1) phases supplied for an
    independent phase-gauge stress test. No velocity or energy integral used.
    """
    H=np.asarray(H)
    if H.ndim!=4 or H.shape[-2:]!=(2,2) or H.shape[0]!=H.shape[1]:
        raise ValueError('a periodic square k grid and two levels required')
    _,u=np.linalg.eigh(H)
    state=u[..., :,band]
    if phases is not None:
        phases=np.asarray(phases)
        if phases.shape!=H.shape[:2]:raise ValueError('phase shape mismatch')
        state=state*np.exp(1j*phases)[...,None]
    ux=np.einsum('...a,...a->...',state.conj(),np.roll(state,-1,axis=0))
    uy=np.einsum('...a,...a->...',state.conj(),np.roll(state,-1,axis=1))
    if min(abs(ux).min(),abs(uy).min())<1e-12:raise ValueError('overlap link vanishes')
    ux/=abs(ux);uy/=abs(uy)
    return float(np.sum(np.angle(ux*np.roll(uy,-1,axis=0)/(
        np.roll(ux,-1,axis=1)*uy)))/(2*np.pi))


def symmetric_spectral_xx(H,vx,mu=0.,eta=.05):
    """Finite-eta clean dissipative Fermi-surface control, NOT the Hall sea."""
    e,u=np.linalg.eigh(H);weight=eta/((mu-e)**2+eta**2)
    U=u@np.einsum('...i,ij->...ij',weight,I)@u.swapaxes(-1,-2).conj()
    return float(np.real(np.mean(np.trace(vx@U@vx@U,axis1=-2,axis2=-1))))
