"""Mixed layer/orbital response of a declared effective Bloch Hamiltonian.

Natural units e=hbar=1. Fields u,b are independent *probe coordinates*:
H=H0+u P-b M+u*b C. All matrices use one fixed orthonormal orbital basis.
Jets store (value, d/du, d/db, d2/du db); slot 3 has no factorial.
No Lorentz force, magnetic phase-space drift correction, or disorder vertex.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import numpy as np

PAIRS = (((0, 0),), ((1, 0), (0, 1)), ((2, 0), (0, 2)),
         ((3, 0), (1, 2), (2, 1), (0, 3)))


def jet(value, a=None, b=None, ab=None):
    value = np.asarray(value, complex)
    return np.stack([value] + [np.zeros_like(value) if v is None else
                              np.broadcast_to(np.asarray(v, complex), value.shape)
                              for v in (a, b, ab)])


def mul(a, b):
    return np.stack([sum(a[i] * b[j] for i, j in pairs) for pairs in PAIRS])


def mm(a, b):
    return np.stack([sum(a[i] @ b[j] for i, j in pairs) for pairs in PAIRS])


def adj(a):
    return a.conj().swapaxes(-1, -2)


def reciprocal(a):
    v, x, y, xy = a
    return np.stack((1/v, -x/v**2, -y/v**2, 2*x*y/v**3-xy/v**2))


def _hermitian(a, name):
    if not np.isfinite(a).all() or not np.allclose(a, a.conj().swapaxes(-1, -2),
                                                atol=1e-11, rtol=1e-10):
        raise ValueError(f'{name} must be finite Hermitian matrices')


@dataclass
class Fields:
    """H0,P,M,C and their analytic k derivatives, with k shape (Nk,2).

    h,p,m,c: (Nk,N,N); dh,dp,dm,dc: (2,Nk,N,N). C defaults to zero.
    M is supplied explicitly; it is NOT inferred from layer or band indices.
    """
    h: np.ndarray
    dh: np.ndarray
    p: np.ndarray
    dp: np.ndarray
    m: np.ndarray
    dm: np.ndarray
    c: np.ndarray | None = None
    dc: np.ndarray | None = None

    def arrays(self):
        h = np.asarray(self.h, complex)
        if h.ndim != 3 or h.shape[-1] != h.shape[-2] or h.shape[-1] < 2:
            raise ValueError('h must have shape (Nk,N,N), N>=2')
        if h.shape[0] == 0:
            raise ValueError('empty k sample')
        if (self.c is None) != (self.dc is None):
            raise ValueError('supply both C and dC, including explicit zero dC for constant C')
        values = [h, self.dh, self.p, self.dp, self.m, self.dm,
                  np.zeros_like(h) if self.c is None else self.c,
                  np.zeros((2,)+h.shape, complex) if self.dc is None else self.dc]
        out = []
        for index, value in enumerate(values):
            a = np.asarray(value, complex)
            shape = h.shape if index % 2 == 0 else (2,)+h.shape
            if a.shape != shape:
                raise ValueError(f'field array {index}: expected {shape}, got {a.shape}')
            _hermitian(a, f'field array {index}')
            out.append(a)
        return out

    def jets(self):
        h, dh, p, dp, m, dm, c, dc = self.arrays()
        return jet(h, p, -m, c), jet(dh, dp, -dm, dc)

    def at(self, u, b):
        if not np.isfinite([u, b]).all():
            raise ValueError('finite field coordinates required')
        h, dh, p, dp, m, dm, c, dc = self.arrays()
        return h+u*p-b*m+u*b*c, dh+u*dp-b*dm+u*b*dc


def eigenjets(hj, gap_tol=1e-8):
    """Normalized nondegenerate Rayleigh--Schrodinger derivatives, both orders.

    Includes |u_ub> on BOTH sides, first-order cross terms and normalization.
    Rejects degeneracies rather than using a gauge-dependent single-band limit.
    """
    hj = np.asarray(hj, complex)
    if hj.ndim != 4 or hj.shape[0] != 4 or hj.shape[-1] != hj.shape[-2]:
        raise ValueError('H jet must be (4,Nk,N,N)')
    if hj.shape[1] == 0 or hj.shape[-1] < 2:
        raise ValueError('nonempty k sample and at least two bands required')
    _hermitian(hj, 'H jet')
    if not np.isfinite(gap_tol) or gap_tol <= 0:
        raise ValueError('positive gap_tol required')
    e, u = np.linalg.eigh(hj[0])
    gap = float(np.min(np.diff(e, axis=-1)))
    if gap <= gap_tol:
        raise ValueError('unresolved band degeneracy: use a separate projector protocol')
    ud = u.conj().swapaxes(-1, -2)
    a, b, c = (ud @ hj[i] @ u for i in (1, 2, 3))
    diagonal = np.eye(e.shape[-1], dtype=bool)
    # Row m, column n: epsilon_n-epsilon_m.
    d = e[:, None, :]-e[:, :, None]
    d = np.where(diagonal, 1., d)
    ua = np.where(diagonal, 0, a/d)
    ub = np.where(diagonal, 0, b/d)
    ea = np.diagonal(a, axis1=-2, axis2=-1).real
    eb = np.diagonal(b, axis1=-2, axis2=-1).real
    source = c+a@ub+b@ua
    eab = np.diagonal(source, axis1=-2, axis2=-1).real
    uab = (source-ua*eb[:, None, :]-ub*ea[:, None, :])/d
    norm = -np.real(np.sum(ua.conj()*ub, axis=-2))
    uab = np.where(diagonal, norm[:, None, :], uab)
    uj = jet(u, u@ua, u@ub, u@uab)
    return jet(e, ea, eb, eab), uj, gap


def geometry(hj, vj, gap_tol=1e-8):
    """Mixed energy, Berry curvature, Berry connection and group velocity."""
    hj, vj = np.asarray(hj, complex), np.asarray(vj, complex)
    if vj.shape != (4, 2)+hj.shape[1:]:
        raise ValueError('velocity jet must be (4,2,Nk,N,N)')
    _hermitian(vj, 'velocity jet')
    ej, uj, gap = eigenjets(hj, gap_tol)
    vb = np.stack([mm(mm(adj(uj), vj[:, i]), uj) for i in range(2)], axis=1)
    d = ej[:, :, None, :]-ej[:, :, :, None]
    diagonal = np.eye(d.shape[-1], dtype=bool)
    d[:, :, diagonal] = np.array([1, 0, 0, 0])[:, None, None]
    inv = reciprocal(d)
    r = np.stack([1j*mul(vb[:, i], inv) for i in range(2)], axis=1)
    r[:, :, :, diagonal] = 0
    omega = -2*np.imag(mul(r[:, 0], r[:, 1].swapaxes(-1, -2)).sum(axis=-1))
    velocity = np.diagonal(vb, axis1=-2, axis2=-1).real
    return dict(energy=ej.real, omega=omega, velocity=velocity, r=r,
                eigenvectors=uj, band_velocity=vb, min_direct_gap=gap)


def orbital_moment(h, dh, ddh=None, gap_tol=1e-8):
    """Eq. A5 operator in orbital basis, including diagonal A1 continuation.

    With ddh[i,j]=d_i d_j H, return analytic d_k M too. This is a REFERENCE
    operator: freezing it in Fields specifies a model, not full magnetic
    minimal coupling or a gate-recomputed self-consistent wave-packet moment.
    """
    h, dh = np.asarray(h, complex), np.asarray(dh, complex)
    if ddh is None:
        ddh = np.zeros((2, 2)+h.shape, complex)
        derivative_requested = False
    else:
        ddh = np.asarray(ddh, complex)
        if ddh.shape != (2, 2)+h.shape:
            raise ValueError('ddh must be (2,2,Nk,N,N)')
        _hermitian(ddh, 'H Hessian')
        if not np.allclose(ddh[0, 1], ddh[1, 0], rtol=1e-9, atol=1e-11):
            raise ValueError('H Hessian must have commuting partial derivatives')
        derivative_requested = True
    g = geometry(jet(h, dh[0], dh[1]), jet(dh, ddh[:, 0], ddh[:, 1]), gap_tol)
    e, u, v = g['energy'], g['eigenvectors'], g['band_velocity']
    nband = h.shape[-1]
    mb = np.zeros_like(u)
    for n in range(nband):
        for p in range(nband):
            for l in range(nband):
                if l in (n, p):
                    continue
                denominator = reciprocal(e[:, :, l]-e[:, :, n])+reciprocal(e[:, :, l]-e[:, :, p])
                cross = (mul(v[:, 0, :, n, l], v[:, 1, :, l, p])
                         - mul(v[:, 1, :, n, l], v[:, 0, :, l, p]))
                mb[:, :, n, p] += -.25j*mul(denominator, cross)
    mo = mm(mm(u, mb), adj(u))
    _hermitian(mo[:3], 'orbital moment and first derivatives')
    return mo[0], (mo[1:3] if derivative_requested else None)


def fermi_derivatives(e, mu, temperature):
    """Energy derivatives f,f',f'',f''' with kBT in MODEL ENERGY units."""
    if not np.isfinite(mu) or not np.isfinite(temperature) or temperature < 0:
        raise ValueError('finite mu and nonnegative thermal energy required')
    e = np.asarray(e, float)
    if temperature == 0:
        f = (e < mu).astype(float)
        if np.any(np.abs(e-mu) < 1e-9) or np.any(np.ptp(f, axis=0) != 0):
            raise ValueError('T=0 only for globally filled/empty sampled bands; metal needs T>0')
        return f, np.zeros_like(f), np.zeros_like(f), np.zeros_like(f)
    f = 1/(1+np.exp(np.clip((e-mu)/temperature, -700, 700)))
    g = f*(1-f)
    return f, -g/temperature, g*(1-2*f)/temperature**2, -g*(1-6*f+6*f*f)/temperature**3


def _weights(w, nk):
    w = np.asarray(w, float)
    if w.shape != (nk,) or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError('explicit positive finite integration weights of shape (Nk,) required')
    return w


def transport(g, weights, mu, temperature):
    """Eq6 / A16 Berry-current and Eq7 / full BAND-DRIFT decomposition.

    Full means the mixed derivative of these two stated functionals, not all
    magnetotransport mechanisms. Extrinsic coefficients are per tau.
    """
    e, o, v = g['energy'], g['omega'], g['velocity']
    w = _weights(weights, len(e[0]))
    f, fp, fpp, fppp = fermi_derivatives(e[0], mu, temperature)
    integ = lambda a: float(np.sum(w[:, None]*np.real(a)))
    terms = dict(sea=integ(o[3]*f), energy_surface=integ(o[0]*e[3]*fp),
                 electric_geometry_magnetic_energy=integ(o[1]*e[2]*fp),
                 magnetic_geometry_electric_energy=integ(o[2]*e[1]*fp),
                 occupation_cross=integ(o[0]*e[1]*e[2]*fpp))
    eq6 = terms['sea']+terms['energy_surface']
    full = sum(terms.values())
    eq7, drift = np.zeros((2, 2)), np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            va, vb = v[:, a], v[:, b]
            eq7[a, b] = -integ((va[3]*vb[0]+va[0]*vb[3])*fp+va[0]*vb[0]*e[3]*fpp)
            fj = jet(fp, e[1]*fpp, e[2]*fpp, e[3]*fpp+e[1]*e[2]*fppp)
            drift[a, b] = -integ(mul(mul(va, vb), fj)[3])
    hall_tensor = lambda h: np.array([[0., -h], [h, 0.]])
    return dict(intrinsic_eq6=eq6, intrinsic_full=full, intrinsic_terms=terms,
                intrinsic_tensor=hall_tensor(full), extrinsic_eq7=eq7,
                band_drift_full=drift, band_drift_cross=drift-eq7,
                min_direct_gap=g['min_direct_gap'])


def exact_geometry(h, dh, gap_tol=1e-8):
    """Non-jet spectral oracle, for independent finite-field differencing."""
    h, dh = np.asarray(h, complex), np.asarray(dh, complex)
    if h.ndim != 3 or h.shape[-1] != h.shape[-2] or dh.shape != (2,)+h.shape:
        raise ValueError('oracle requires H(Nk,N,N) and dH(2,Nk,N,N)')
    _hermitian(h, 'H'); _hermitian(dh, 'dH')
    e, u = np.linalg.eigh(h)
    if np.min(np.diff(e, axis=-1)) <= gap_tol:
        raise ValueError('unresolved band degeneracy in finite-field oracle')
    v = np.stack([u.conj().swapaxes(-1, -2) @ dh[i] @ u for i in range(2)])
    delta = e[:, :, None]-e[:, None, :]
    diagonal = np.eye(e.shape[-1], dtype=bool)
    delta = np.where(diagonal, np.inf, delta)
    o = -2*np.imag(v[0]*v[1].swapaxes(-1, -2)/delta**2).sum(axis=-1)
    return e, o, np.diagonal(v, axis1=-2, axis2=-1).real


def exact_transport(fields, weights, mu, temperature, u, b):
    e, o, v = exact_geometry(*fields.at(u, b))
    w = _weights(weights, len(e))
    f, fp, _, _ = fermi_derivatives(e, mu, temperature)
    hall = float(np.sum(w[:, None]*o*f))
    drift = -np.einsum('k,akn,bkn,kn->ab', w, v, v, fp)
    return hall, drift


def mixed_difference(fn: Callable, hu=1e-3, hb=1e-3):
    if not np.isfinite([hu, hb]).all() or min(hu, hb) <= 0:
        raise ValueError('positive finite field steps required')
    return (np.asarray(fn(hu, hb))-np.asarray(fn(hu, -hb))
            -np.asarray(fn(-hu, hb))+np.asarray(fn(-hu, -hb)))/(4*hu*hb)


def projector_curvature(hamiltonian: Callable, points, step=1e-4):
    """Independent k-difference projector formula iTr P[DxP,DyP]."""
    if not np.isfinite(step) or step <= 0:
        raise ValueError('positive step required')
    def project(k):
        e, u = np.linalg.eigh(hamiltonian(k))
        if np.min(np.diff(e, axis=-1)) <= 1e-8:
            raise ValueError('unresolved degeneracy in projector oracle')
        return np.einsum('kin,kjn->knij', u, u.conj())
    p = project(points)
    dp = []
    for a in range(2):
        offset = np.eye(2)[a]*step
        dp.append((project(points+offset)-project(points-offset))/(2*step))
    return np.real(1j*np.trace(p@(dp[0]@dp[1]-dp[1]@dp[0]), axis1=-2, axis2=-1))
