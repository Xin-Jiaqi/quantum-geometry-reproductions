"""Zou et al. (PRL 135, 226901, 2025), SM Secs II.2 and IV.5.

Dimensionless lattice momentum, energies in the SM tight-binding energy unit,
hbar=1.  This is a *leading DBC singularity* audit, not a full SHG material engine.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple
import numpy as np

C0, C1, C2, T, D1, D2, EG = -0.75, 0.5, 1.0, 1.0, 0.3, 0.3, -6.0
EPSILONS = (0.75, 0.90, 1.03, 1.50, 2.00, 3.00)
SOURCE_E0 = (6.768, 6.919, 7.049, 7.520, 8.021, 9.022)
SOURCE_VRATIO = (0.701, 0.849, 0.978, 1.445, 1.941, 2.935)


def matrices(kx: float, ky: float, epsilon1: float):
    """Return H, (dHx,dHy), (dxxH,dxyH,dyyH) analytically.

    The basis is the SM Eq. (IV.5) (p_y, p_x, p_y); the global phase of
    orbital hopping and the third-band onsite term are kept as printed.
    """
    x, y = float(kx), float(ky)
    a = C0 + C1 * np.cos(y) + C2 * np.cos(x)
    h = np.array([
        [a, -1j * T * np.sin(x), -1j * D1 * np.sin(y)],
        [0, -a, -1j * D2 * np.sin(x)],
        [0, 0, EG + epsilon1 * (C0 + C1 * np.cos(y))],
    ], dtype=complex)
    vx = np.array([
        [-C2 * np.sin(x), -1j * T * np.cos(x), 0],
        [0, C2 * np.sin(x), -1j * D2 * np.cos(x)],
        [0, 0, 0],
    ], dtype=complex)
    vy = np.array([
        [-C1 * np.sin(y), 0, -1j * D1 * np.cos(y)],
        [0, C1 * np.sin(y), 0],
        [0, 0, -epsilon1 * C1 * np.sin(y)],
    ], dtype=complex)
    wxx = np.array([
        [-C2 * np.cos(x), +1j * T * np.sin(x), 0],
        [0, C2 * np.cos(x), +1j * D2 * np.sin(x)],
        [0, 0, 0],
    ], dtype=complex)
    wyy = np.array([
        [-C1 * np.cos(y), 0, +1j * D1 * np.sin(y)],
        [0, C1 * np.cos(y), 0],
        [0, 0, -epsilon1 * C1 * np.cos(y)],
    ], dtype=complex)
    def hermitian(a):
        return a + np.triu(a, 1).conj().T
    return hermitian(h), (hermitian(vx), hermitian(vy)), (hermitian(wxx), np.zeros((3, 3), complex), hermitian(wyy))


def nodal_parameters(epsilon1: float):
    """Positive-ky crossing between the top two TB bands; SM Table S1."""
    def signed(y):
        h = matrices(0., y, epsilon1)[0]
        aa, ev = h[0, 0].real, h[2, 2].real
        top = (aa + ev + np.sqrt((aa - ev) ** 2 + 4 * abs(h[0, 2]) ** 2)) / 2
        return float(top + aa)
    lo, hi = 1.3, 2.7
    if signed(lo) * signed(hi) >= 0:
        raise ValueError('Expected mirror nodal point not bracketed')
    for _ in range(70):
        mid = (lo + hi) / 2
        if signed(mid) * signed(lo) <= 0: hi = mid
        else: lo = mid
    y0 = (lo + hi) / 2
    E = np.linalg.eigvalsh(matrices(0, y0, epsilon1)[0])
    def bands(y):
        h = matrices(0, y, epsilon1)[0]
        aa, ev = h[0, 0].real, h[2, 2].real
        gap = np.sqrt((aa - ev) ** 2 + 4 * abs(h[0, 2]) ** 2)
        return np.array(((aa + ev - gap) / 2, (aa + ev + gap) / 2, -aa))
    delta = 1e-5
    v = (bands(y0 + delta) - bands(y0 - delta)) / (2 * delta)
    vc = abs((v[1] - v[2]) / 2)
    vv = abs((v[1] + v[2]) / 2 - v[0])
    return {'epsilon1': float(epsilon1), 'ky': float(y0),
            'epsilon0': float(E[1] - E[0]), 'vc': float(vc),
            'vv': float(vv), 'vv_over_vc': float(vv / vc)}


class BandData(NamedTuple):
    energies: np.ndarray
    berry: np.ndarray  # r[n,m,cart], r_nm=-i<v_n|dH|v_m>/(E_n-E_m)
    velocity: np.ndarray  # v[cart,n,m]
    w: np.ndarray  # w[cart,cart,n,m]


def band_data(kx: float, ky: float, epsilon1: float) -> BandData:
    h, v, w = matrices(kx, ky, epsilon1)
    E, U = np.linalg.eigh(h)
    V = np.stack([U.conj().T @ x @ U for x in v])
    W = np.array([[U.conj().T @ w[0] @ U, U.conj().T @ w[1] @ U],
                  [U.conj().T @ w[1] @ U, U.conj().T @ w[2] @ U]])
    R = np.zeros((3, 3, 2), complex)
    for n in range(3):
        for m in range(3):
            if n != m:
                if abs(E[n] - E[m]) < 1e-13:
                    raise ValueError('Abelian Berry connection undefined at exact degeneracy')
                R[n, m] = -1j * V[:, n, m] / (E[n] - E[m])
    return BandData(E, R, V, W)


def covariant_derivative(d: BandData, n: int, m: int, deriv: int, optical: int):
    """Full SM Eq. (S10) in hbar=1 units, for n!=m and nondegenerate k."""
    if n == m: raise ValueError('Requires n != m')
    E, V, W = d.energies, d.velocity, d.w
    om = E[n] - E[m]
    delta_d = V[deriv, n, n] - V[deriv, m, m]
    delta_o = V[optical, n, n] - V[optical, m, m]
    q = (V[optical, n, m] * delta_d + V[deriv, n, m] * delta_o) / om - W[deriv, optical, n, m]
    for p in range(3):
        if p not in (n, m):
            q += V[optical, n, p] * V[deriv, p, m] / (E[p] - E[m])
            q -= V[deriv, n, p] * V[optical, p, m] / (E[n] - E[p])
    return 1j * q / om


def hermitian_connection(d: BandData, n: int, m: int, a: int, b: int, c: int):
    """C_{nm}^{abc} = r_nm^b D^a[r_mn^c] (note the reversed indices)."""
    return d.berry[n, m, b] * covariant_derivative(d, m, n, a, c)


def singular_connection(r: np.ndarray, n: int, m: int, a: int, b: int, c: int):
    """Leading term of C_nm^{abc} from SM (S12), for c,c'={1,2}, v=0.

    Note the derivative acts on r_{mn} while its prefactor is r_{nm}.
    Reversing this order produces a spurious non-cancellation.
    """
    z = 0j
    if m in (1, 2):
        partner = 3 - m
        z += 1j * r[n, m, b] * r[m, partner, a] * r[partner, n, c]
    if n in (1, 2):
        partner = 3 - n
        z -= 1j * r[n, m, b] * r[m, partner, c] * r[partner, n, a]
    return z


@dataclass(frozen=True)
class SHGPieces:
    inter_2omega: complex
    inter_omega_a: complex
    inter_omega_b: complex
    mixed_2omega: complex
    mixed_omega_a: complex
    mixed_omega_b: complex

    @property
    def interband(self):
        return self.inter_2omega + self.inter_omega_a + self.inter_omega_b

    @property
    def mixed_singular(self):
        return self.mixed_2omega + self.mixed_omega_a + self.mixed_omega_b

    @property
    def residual(self):
        return self.interband + self.mixed_singular

    @property
    def double_residual(self):
        return self.inter_2omega + self.mixed_2omega

    @property
    def single_residual(self):
        return self.inter_omega_a + self.inter_omega_b + self.mixed_omega_a + self.mixed_omega_b


def shg_pieces(energies: np.ndarray, berry: np.ndarray, frequency: float,
               axes=(0, 0, 1)) -> SHGPieces:
    """Direct SM (S6)+(S19); omit common e^3/(2 hbar^2), k-space integration.

    SM (S19) includes only the leading divergent Hermitian-connection terms.
    Consequently `residual` is finite as the DBC is approached, but is *not*
    the full finite SHG integrand. Frequency must be off an optical pole.
    """
    a, b, g = axes
    e, r = energies, berry
    f = np.array([1., 0., 0.])
    e2 = e1a = e1b = 0j
    for n in range(3):
        for m in range(3):
            for l in range(3):
                if len({n, m, l}) != 3: continue
                wln, wml, wmn = e[l] - e[n], e[m] - e[l], e[m] - e[n]
                triple = r[n, m, a] * (r[m, l, b] * r[l, n, g] + r[m, l, g] * r[l, n, b]) / (wln - wml)
                e2 += triple * 2 * (f[n] - f[m]) / (wmn - 2 * frequency)
                e1a += triple * (f[l] - f[n]) / (wln - frequency)
                e1b += triple * (f[m] - f[l]) / (wml - frequency)
    i2 = i1a = i1b = 0j
    for n in range(3):
        for m in range(3):
            if n == m or f[n] == f[m]: continue
            wmn = e[m] - e[n]
            sc = lambda p, q, t, u, v: singular_connection(r, p, q, t, u, v)
            A = sc(n, m, g, a, b) + sc(n, m, b, a, g)
            B = sc(m, n, g, b, a) + sc(m, n, b, g, a)
            C = sc(m, n, a, g, b) + sc(m, n, a, b, g)
            i2 += 1j * (f[n] - f[m]) * 2 * A / (wmn * (wmn - 2 * frequency))
            i1a += 1j * (f[n] - f[m]) * B / (wmn * (wmn - frequency))
            i1b -= 1j * (f[n] - f[m]) * C / (2 * wmn * (wmn - frequency))
    return SHGPieces(e2, e1a, e1b, i2, i1a, i1b)


def benchmark(epsilon1: float, q: float, axes=(0, 0, 1), frequency_fraction=0.7):
    """Generic oblique approach to node; avoid a high-symmetry zero direction."""
    info = nodal_parameters(epsilon1)
    d = band_data(q, info['ky'] + 0.37 * q, epsilon1)
    p = shg_pieces(d.energies, d.berry, frequency_fraction * info['epsilon0'], axes)
    return {'epsilon1': float(epsilon1), 'q': float(q), 'tensor': ''.join('xy'[i] for i in axes),
            'photon_energy': frequency_fraction * info['epsilon0'],
            'inter_abs': float(abs(p.interband)), 'mixed_abs': float(abs(p.mixed_singular)),
            'residual_abs': float(abs(p.residual)),
            'double_residual_abs': float(abs(p.double_residual)),
            'single_residual_abs': float(abs(p.single_residual)),
            'residual_over_inter': float(abs(p.residual) / abs(p.interband))}

