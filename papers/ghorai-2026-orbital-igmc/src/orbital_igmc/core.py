"""Low-frequency orbital magnetic displacement current, with source audit.

Units: q=k*a dimensionless, H and dH/dq in eV, b=B/B_unit,
B_unit=hbar/(e*a**2). M is the energy coupling to b in eV.
Every coefficient uses B(t)=B0*cos(omega*t) and the +omega Fourier
coefficient divided by B0 (not the conventional full phasor).

The two REQUIRED protocols are kept distinct: printed_D10 implements
all three orbital terms of the uploaded PRB; density_A9 evaluates its
kinetic equation directly. They do not numerically agree; see SOURCE_AUDIT.md.
No collision, resonant, Fermi-surface transport, or Maxwell circuit solver.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

E_CHARGE = 1.602176634e-19
HBAR_EV_S = 6.582119569e-16
KB_EV_K = 8.617333262145e-5
MU_B_EV_T = 5.7883818060e-5
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1, -1]).astype(complex)
TX, TZ = np.kron(SX, I2), np.kron(SZ, I2)
IX, IY, IZ = (np.kron(I2, s) for s in (SX, SY, SZ))
I4 = np.eye(4, dtype=complex)


@dataclass(frozen=True)
class CuMnAs:
    """PRB Eq.(37). The missing in-plane length MUST be supplied explicitly."""
    a_angstrom: float
    phi: float = np.pi / 12
    t0: float = 0.08
    t1: float = 1.0
    alpha: float = 0.8
    exchange: float = 0.85
    g: float = 2.0
    group_sizes: tuple = (2, 2)

    def __post_init__(self):
        vals = [self.a_angstrom, self.phi, self.t0, self.t1,
                self.alpha, self.exchange, self.g]
        if not np.isfinite(vals).all() or self.a_angstrom <= 0 or self.g < 0:
            raise ValueError('Finite parameters, a_angstrom>0 and g>=0 required')
        if self.group_sizes != (2, 2):
            raise ValueError('CuMnAs has two rank-two PT doublets')

    @property
    def field_unit_T(self):
        # hbar[J s]/e[C] = hbar[eV s] numerically.
        return HBAR_EV_S / (self.a_angstrom * 1e-10)**2

    def operators(self, q):
        q = np.asarray(q, float)
        if q.ndim != 2 or q.shape[1] != 2 or not np.isfinite(q).all():
            raise ValueError('q must be finite (Nk,2) dimensionless momenta')
        x, y = q.T
        def lift(s): return s[:, None, None]
        c = -self.t0 * (np.cos(x) + np.cos(y))
        v = -2 * self.t1 * np.cos(x / 2) * np.cos(y / 2)
        hx = self.exchange * np.cos(self.phi) - self.alpha * np.sin(y)
        hy = self.exchange * np.sin(self.phi) + self.alpha * np.sin(x)
        H = lift(c)*I4 + lift(v)*TX + lift(hx)*(TZ@IX) + lift(hy)*(TZ@IY)
        dx = (lift(self.t0*np.sin(x))*I4
              + lift(self.t1*np.sin(x/2)*np.cos(y/2))*TX
              + lift(self.alpha*np.cos(x))*(TZ@IY))
        dy = (lift(self.t0*np.sin(y))*I4
              + lift(self.t1*np.cos(x/2)*np.sin(y/2))*TX
              - lift(self.alpha*np.cos(y))*(TZ@IX))
        M = -0.5*self.g*MU_B_EV_T*self.field_unit_T*IZ
        return H, np.stack([dx, dy]), M


def mesh(n):
    """Uniform midpoint BZ; mean represents integral d2q/(2pi)^2."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 8:
        raise ValueError('Integer mesh size >=8 required')
    axis = -np.pi + (np.arange(n)+0.5)*2*np.pi/n
    return np.stack(np.meshgrid(axis, axis, indexing='ij'), -1).reshape(-1, 2)


@dataclass
class Geometry:
    energies: np.ndarray
    velocities: np.ndarray
    connection: np.ndarray
    spin_moment: np.ndarray
    inverse_gap: np.ndarray
    metric: np.ndarray
    orbital_moment_term: np.ndarray
    orbital_multiband_term: np.ndarray
    spin_term: np.ndarray
    group_sizes: tuple
    min_gap_eV: float


def geometry(H, derivatives, spin_moment, *, group_sizes, eigenvectors=None,
             degeneracy_tol=1e-9, gap_floor=1e-6):
    """Spectral blocks for a fixed orthonormal orbital basis (B along z).

    group_sizes is explicit, e.g. (2,2) or (1,1,1). Degenerate groups must
    persist locally: projected dH/dq must be scalar within each group.
    Near crossings and non-scalar intragroup derivatives fail closed.
    Individual degenerate-band entries are gauge dependent; only occupied
    block traces are returned by integrate. No arbitrary gap regularization.
    """
    H, V, M = map(lambda a: np.asarray(a, complex), (H, derivatives, spin_moment))
    if H.ndim != 3 or H.shape[-1] != H.shape[-2] or H.shape[0] == 0:
        raise ValueError('H must have shape (Nk,Nb,Nb)')
    nk, nb, _ = H.shape
    sizes = tuple(group_sizes)
    if (len(sizes) < 2 or any(isinstance(s, bool) or not isinstance(s, (int,np.integer))
                             or s < 1 for s in sizes) or sum(sizes) != nb):
        raise ValueError('group_sizes must partition all bands into >=2 groups')
    if V.shape != (2, nk, nb, nb) or M.shape not in [(nb,nb), H.shape]:
        raise ValueError('dH shape=(2,Nk,Nb,Nb); M shape=(Nb,Nb) or H.shape')
    if not (0 < degeneracy_tol < gap_floor) or not np.isfinite([degeneracy_tol,gap_floor]).all():
        raise ValueError('Require 0<degeneracy_tol<gap_floor')
    for op in (H,V,M):
        if not np.isfinite(op).all() or not np.allclose(op,op.swapaxes(-1,-2).conj(),rtol=0,atol=2e-11):
            raise ValueError('Finite Hermitian operators required')
    E, U = np.linalg.eigh(H)
    if eigenvectors is not None:
        U = np.asarray(eigenvectors, complex)
        if (U.shape != H.shape or not np.isfinite(U).all()
                or not np.allclose(U.swapaxes(-1,-2).conj()@U,np.eye(nb),atol=2e-10,rtol=0)
                or not np.allclose(H@U,U*E[:,None,:],atol=2e-9,rtol=0)):
            raise ValueError('Eigenvectors must be orthonormal and match sorted energies')
    Ud = U.swapaxes(-1,-2).conj()
    V = Ud@V@U
    M = Ud@M@U
    start = 0
    for size in sizes:
        sl = slice(start,start+size)
        if np.max(np.ptp(E[:,sl],axis=-1)) > degeneracy_tol:
            raise ValueError('Declared degenerate group is split')
        block = V[:,:,sl,sl]
        mean = np.trace(block,axis1=-2,axis2=-1)/size
        if not np.allclose(block,mean[...,None,None]*np.eye(size),atol=2e-9,rtol=0):
            raise ValueError('Degeneracy splits locally; scalar intragroup derivatives required')
        start += size
    labels = np.repeat(np.arange(len(sizes)),sizes)
    between = labels[:,None] != labels[None,:]
    delta = E[:,None,:]-E[:,:,None]  # E_m-E_n
    gap = float(np.min(np.abs(delta[:,between])))
    if gap < gap_floor:
        raise ValueError('Unresolved intergroup gap; no intrinsic denominator at this crossing')
    inv = np.divide(1.,delta,out=np.zeros_like(delta),where=between[None])
    R = 1j*V*inv
    vd = np.diagonal(V,axis1=-2,axis2=-1).real
    vx, vy = V; rx, ry = R
    # M_mn=sum_{p outside n group} M_mpn, from PRB Eqs.(11),(17).
    ML = -.5*(vx@ry-vy@rx + ry*vd[0,:,None,:]-rx*vd[1,:,None,:])
    eye = np.eye(nb)[None,:,:,None]
    # Axes: k,m,n,p. R is zero WITHIN a spectral group.
    ML3 = -.5*((vx[:,:,:,None]+eye*vd[0,:,None,None,:])*ry[:,None,:,:]
              -(vy[:,:,:,None]+eye*vd[1,:,None,None,:])*rx[:,None,:,:])
    moment = -np.real(np.sum(ML.swapaxes(-1,-2)[None]*R*inv,axis=-1))
    multi = np.array([np.einsum('kmnp,kpm,kmp->kn',ML3,ra,-inv).real for ra in R])
    spin = np.real(np.sum(M.swapaxes(-1,-2)[None]*R*inv,axis=-1))
    metric = np.array([[np.sum((rc*ra.swapaxes(-1,-2)).real,axis=-1)
                        for ra in R] for rc in R])
    return Geometry(E,V,R,M,inv,metric,moment,multi,spin,sizes,gap)


def _occupation(E, mu, temperature_K):
    if not np.isfinite([mu,temperature_K]).all() or temperature_K < 0:
        raise ValueError('Finite mu and temperature_K>=0 required')
    if temperature_K == 0:
        # The metric derivative in a T=0 metal needs a Fermi-line method.
        if np.any((E.min(axis=0) <= mu) & (mu <= E.max(axis=0))):
            raise ValueError('Zero-T metal not supported; use finite T or a Fermi-line solver')
        return (E<mu).astype(float), np.zeros_like(E)
    x = (E-mu)/(KB_EV_K*temperature_K)
    tail = np.exp(-np.abs(x))
    f = np.where(x>=0,tail/(1+tail),1/(1+tail))
    fp = -tail/((1+tail)**2*KB_EV_K*temperature_K)
    return f, fp


def integrate(geom, mu_eV, temperature_K, *, protocol):
    """BZ-averaged positive-harmonic polarization kernel P_+/B0=-e G/(a Bunit).

    Call on a complete uniform periodic BZ. Printed metric derivative is
    integrated by parts: -int f curl(g)/2 = int grad(f) cross(g)/2.
    density_A9 is evaluated from its driving matrix, NOT by rescaling D10.
    """
    if protocol not in ('printed_D10','density_A9'):
        raise ValueError('Explicit protocol required: printed_D10 or density_A9')
    E, V, R, inv = geom.energies, geom.velocities, geom.connection, geom.inverse_gap
    f, fp = _occupation(E,mu_eV,temperature_K)
    vd = np.diagonal(V,axis1=-2,axis2=-1).real
    def bz(a): return np.mean(np.sum(a,axis=-1),axis=-1)
    if protocol == 'printed_D10':
        # epsilon_zbc with (b,c)=(x,y),(y,x).
        metric = bz(.5*fp[None]*(vd[0][None]*geom.metric[1]-vd[1][None]*geom.metric[0]))
        moment = bz(f[None]*geom.orbital_moment_term)
        multi = bz(f[None]*geom.orbital_multiband_term)
        spin = bz(f[None]*geom.spin_term)
        orbital = metric+moment+multi
        parts = dict(metric=metric,moment=moment,multiband=multi)
    else:
        # Frechet derivative d f(H), including equal-energy f' blocks.
        fnm = f[:,:,None]-f[:,None,:]
        divided = -fnm*inv
        same = inv == 0
        divided = np.where(same,fp[:,:,None],divided)
        dx,dy = V*divided
        vx,vy = V
        drive_L = (vy@dx+dx@vy-vx@dy-dy@vx)/4
        drive_S = -.5j*geom.spin_moment*fnm
        def coefficient(drive):
            rho = 1j*drive*inv  # 1/[i(E_n-E_m)], off-group only.
            return np.array([np.mean(np.einsum('kij,kji->k',rho,r).real) for r in R])
        orbital,spin = coefficient(drive_L),coefficient(drive_S)
        parts = {}
    return dict(protocol=protocol,orbital=orbital,spin=spin,total=orbital+spin,
                min_gap_eV=geom.min_gap_eV,**parts)


def response_si(result, model, drive_energy_eV, *, convention):
    """Return complex sheet current/B in A/(m T); omega=drive_energy_eV/hbar.

    harmonic: j(t)=chi*B0*exp(+iwt)+c.c.; phasor: j(t)=Re[chi*B0*exp(+iwt)].
    No effective 3D thickness is inferred. Validity additionally requires
    omega*tau<<1 and hbar/tau<<the intergroup gap; tau is not a fit parameter.
    """
    if convention not in ('harmonic','phasor'):
        raise ValueError('Explicit harmonic or phasor convention required')
    if not np.isfinite(drive_energy_eV) or abs(drive_energy_eV) > .05*result['min_gap_eV']:
        raise ValueError('Low-frequency expansion requires |hbar*omega|<=0.05*min_gap')
    pref = -1j*E_CHARGE*(drive_energy_eV/HBAR_EV_S)/(model.a_angstrom*1e-10*model.field_unit_T)
    if convention == 'phasor': pref *= 2
    return {k:pref*np.asarray(result[k]) for k in ('orbital','spin','total')}


def density_reference(model, q, mu_eV, temperature_K, *, step=2e-5):
    """Independent fixed-orbital finite differences of H, f(H), and projectors.

    Numerically differentiates projectors, not eigenvector phases; then solves
    Eq.(A9) between spectral groups. Used as a reference, not the fast API.
    """
    if not np.isfinite(step) or step<=0:
        raise ValueError('Positive finite difference step required')
    q = np.asarray(q,float)
    H,_,M = model.operators(q)
    E,U = np.linalg.eigh(H); Ud = U.swapaxes(-1,-2).conj()
    edges = np.cumsum((0,)+model.group_sizes)
    def state(points):
        h,_,_ = model.operators(points)
        e,u = np.linalg.eigh(h); ud=u.swapaxes(-1,-2).conj()
        f,_ = _occupation(e,mu_eV,temperature_K)
        rho=(u*f[:,None,:])@ud
        ps=[u[:,:,edges[s]:edges[s+1]]@ud[:,edges[s]:edges[s+1],:]
            for s in range(len(model.group_sizes))]
        return h,rho,ps
    _,rho0,ps = state(q)
    dh=[]; dr=[]; connections=[]
    for a in (0,1):
        dq=np.zeros_like(q);dq[:,a]=step
        hp,rp,pp=state(q+dq);hm,rm,pm=state(q-dq)
        dh.append((hp-hm)/(2*step));dr.append((rp-rm)/(2*step))
        connections.append(1j*sum((pp[s]-pm[s])@ps[s]/(2*step) for s in range(len(ps))))
    vx,vy=dh;dx,dy=dr
    source_L=(vy@dx+dx@vy-vx@dy-dy@vx)/4
    source_S=.5j*(M@rho0-rho0@M)
    labels=np.repeat(np.arange(len(model.group_sizes)),model.group_sizes)
    mask=(labels[:,None]!=labels[None,:])[None]
    delta=E[:,:,None]-E[:,None,:]
    if np.min(np.abs(delta[:,mask[0]]))<1e-6:
        raise ValueError('Reference encountered unresolved intergroup gap')
    result={}
    for key,source in (('orbital',source_L),('spin',source_S)):
        rhob=np.divide(Ud@source@U,1j*delta,out=np.zeros_like(source),where=mask)
        rhor=U@rhob@Ud
        result[key]=np.array([np.mean(np.einsum('kij,kji->k',rhor,r).real) for r in connections])
    return result
