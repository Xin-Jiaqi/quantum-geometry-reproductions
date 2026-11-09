"""Magnetic injection current of the Ezawa 2D Rashba altermagnet.

The physically *distinct* symmetric optical sector is contracted with the
interband quantum metric g_ab = Re[r^a_{-+} r^b_{+-}]. This is NOT the
antisymmetric Berry-commutator kernel used for circular injection.

Units: hbar=e=tau=1 and W*k0=1, momentum in k0, energy in epsilon0; coefficient
C = sigma_inject / [tau*e**3/(hbar**2*W*k0)] for e>0. No SI output.

For a single positive-energy transition (lower->upper), Eq. (7) implies
C^{c;ab}(E) = -(1/pi) int d^2k (f_- - f_+) (partial_c d)
                            g_ab delta(2*d-E),
where H = a*k^2 I + h(k).sigma and d=|h|.
A single monotonically increasing resonance contour can be integrated exactly
using delta(2d-E). The finite-broadening grid is a separate numerical check.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Params:
    rashba: float = 1.0
    altermagnet: float = 0.10
    gap: float = 0.10
    parabola: float = 0.25  # hbar^2/(2M) in epsilon0/k0^2
    mu: float = 0.0

    def validate(self) -> None:
        if self.rashba < 0:
            raise ValueError("rashba must be nonnegative")
        if self.parabola <= abs(self.altermagnet):
            raise ValueError("use parabola>|altermagnet| for positive large-k lower band")


def band_observables(kx, ky, p: Params):
    """Vectorized energies, interband velocity, and two-band quantum metric.

    All arrays have broadcasted momentum shape, except metric (...,2,2)
    and velocity_difference (...,2). H's scalar parabolic term affects
    occupations only: it cancels from interband matrix elements and Delta-v.
    """
    p.validate()
    x, y = np.broadcast_arrays(np.asarray(kx, float), np.asarray(ky, float))
    lam, J, B = p.rashba, p.altermagnet, p.gap
    h = np.stack((-lam*y, lam*x+J*(x*x-y*y), np.broadcast_to(B, x.shape)), axis=-1)
    dhx = np.stack((np.zeros_like(x), lam+2*J*x, np.zeros_like(x)), axis=-1)
    dhy = np.stack((-lam*np.ones_like(x), -2*J*y, np.zeros_like(x)), axis=-1)
    dh = np.stack((dhx, dhy), axis=-2)
    d2 = np.einsum('...i,...i->...', h, h)
    if np.any(d2 <= 0.0):
        raise ValueError("gapless degeneracy: nondegenerate quantum metric undefined")
    d = np.sqrt(d2)
    hd = np.einsum('...i,...ai->...a', h, dh)
    grad = hd / d[..., None]
    metric = 0.25 * (
        np.einsum('...ai,...bi->...ab', dh, dh)/d2[..., None, None]
        - hd[..., :, None]*hd[..., None, :]/(d2[..., None, None]**2)
    )
    h0 = p.parabola*(x*x+y*y)
    return {
        'lower': h0-d, 'upper': h0+d, 'gap': 2*d,
        'velocity_difference': 2*grad, 'metric': metric,
        'spin_field': h,
    }


def _occupancy(obs, mu):
    """T=0 transition weight f_lower-f_upper (may be 0 or 1 here)."""
    return ((obs['lower'] < mu).astype(float)
            - (obs['upper'] < mu).astype(float))


def contour_injection(photon_energy: float, p: Params = Params(),
                      current: str = 'x', optical: str = 'xx',
                      ntheta: int = 720, kmax: float = 3.0) -> float:
    """Sharp-delta 2D resonant-contour injection; single-root protocol.

    Midpoint quadrature in angle plus vectorized bisection in radius. Rejects
    missing roots, nonmonotonic radial dispersion or a gapless model. For
    more general dispersions or multiple roots use the broadened BZ grid.
    """
    if current not in ('x', 'y') or len(optical) != 2 or any(c not in 'xy' for c in optical):
        raise ValueError("current must be x/y and optical must be two x/y letters")
    if ntheta < 16 or kmax <= 0:
        raise ValueError("ntheta>=16 and positive kmax required")
    if p.rashba == 0:
        raise ValueError("Rashba=0: use grid_injection for the inversion-symmetric limit")
    if abs(p.gap) == 0:
        raise ValueError("B=0 has a degeneracy at Gamma; use a finite gap regulator")
    if photon_energy <= 2*abs(p.gap):
        return 0.0
    theta = 2*np.pi*(np.arange(ntheta)+0.5)/ntheta
    co, si = np.cos(theta), np.sin(theta)
    # Precondition: radial gap is monotone for k <= kmax, so each ray has one root.
    # Reject unsafe parameters rather than silently discard additional roots.
    rcheck = np.linspace(0.0, kmax, 81)[1:, None]
    s = band_observables(rcheck*co, rcheck*si, p)
    radial = 0.5*np.sum(s['velocity_difference']*np.stack((co, si), axis=-1), axis=-1)
    if np.any(radial <= 0):
        raise ValueError("non-monotonic radial gap: contour single-root assumption fails")
    if np.any(s['gap'][-1] <= photon_energy):
        raise ValueError("increase kmax: some resonance roots are outside the domain")
    lo, hi = np.zeros(ntheta), np.full(ntheta, kmax)
    for _ in range(60):
        mid = 0.5*(lo+hi)
        gap_mid = band_observables(mid*co, mid*si, p)['gap']
        lo = np.where(gap_mid < photon_energy, mid, lo)
        hi = np.where(gap_mid >= photon_energy, mid, hi)
    r = 0.5*(lo+hi)
    obs = band_observables(r*co, r*si, p)
    grad = 0.5*obs['velocity_difference']
    dradial = grad[..., 0]*co+grad[..., 1]*si
    ci, ai, bi = 'xy'.index(current), 'xy'.index(optical[0]), 'xy'.index(optical[1])
    weight = _occupancy(obs, p.mu)
    return float(-np.mean(r*grad[..., ci]*obs['metric'][..., ai, bi]
                          *weight/np.abs(dradial)))


def grid_injection(photon_energy: float, p: Params = Params(),
                   current: str = 'x', optical: str = 'xx',
                   nk: int = 701, kmax: float = 1.5,
                   gamma: float = .01) -> float:
    """Independent 2D Cartesian-grid Lorentzian evaluation of Eq. (7).

    Controlled kmax/nk/gamma extrapolation is essential: Lorentzian tails
    make this NOT mathematically identical to the sharp-delta calculation.
    """
    if nk < 21 or gamma <= 0 or kmax <= 0:
        raise ValueError("nk>=21, gamma>0, kmax>0 required")
    if current not in ('x', 'y') or len(optical) != 2 or any(c not in 'xy' for c in optical):
        raise ValueError("invalid tensor component")
    k = -kmax+(np.arange(nk)+0.5)*(2*kmax/nk)
    x, y = np.meshgrid(k, k, indexing='ij')
    obs = band_observables(x, y, p)
    ci, ai, bi = 'xy'.index(current), 'xy'.index(optical[0]), 'xy'.index(optical[1])
    delta = gamma/(np.pi*((obs['gap']-photon_energy)**2+gamma**2))
    integrand = (_occupancy(obs, p.mu) * .5*obs['velocity_difference'][..., ci]
                 * obs['metric'][..., ai, bi] * delta)
    return float(-(2*kmax/nk)**2 / np.pi * np.sum(integrand))


def analytic_xxx_eq17(photon_energy: float, p: Params = Params()) -> float:
    """Printed PRB Eq. (17), O(J); only meaningful for E>2|B| below mu cutoff."""
    if photon_energy <= 0:
        raise ValueError("positive photon energy required")
    b = p.gap / photon_energy
    return -p.altermagnet/(8*p.rashba)*(1+8*b*b-48*b**4)


def analytic_xyy_printed_eq19(photon_energy: float, p: Params = Params()) -> float:
    """Literal paper Eq. (19), retained for an explicit discrepancy audit."""
    if photon_energy <= 0:
        raise ValueError("positive photon energy required")
    b = p.gap / photon_energy
    return -p.altermagnet/(8*p.rashba)*(1-2*b*b+4*b**4)
