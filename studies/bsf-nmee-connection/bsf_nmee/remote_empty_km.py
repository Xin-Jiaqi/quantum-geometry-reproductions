"""KM bilayer plus one spinful empty orbital per layer.

The original KM orbitals, spin operators and hopping displacement convention
are retained. The new orbital lies on A or B and couples across the interface
with the same real-space hopping envelope. No MoSSe data are fitted.
"""
import numpy as np
from nmee_lab.clean_core.core import MU_B

from .bilayer_audit import AuditBilayer, AuditSpec
from .bilayer_protocol import ab_ba_shift
from .mixing_channels import adj, split_channels


def ingredients(k, mass=.045, fraction=.25, site=0, angle=np.pi/3, shift=None, rashba=.020):
    unit = AuditBilayer(AuditSpec(shift=ab_ba_shift(fraction) if shift is None else shift,
                                 staggered_eV=mass, tunneling_eV=1., rashba_soc_eV=rashba))
    direction = np.array([np.cos(angle), np.sin(angle)])
    old = unit.matrices(k)
    vold = sum(direction[j] * unit.matrices(k, axis=j) for j in range(2))
    h, v, low, dlow, extra, dextra = [np.zeros((len(k), 12, 12), complex) for _ in range(6)]
    for offset, oldoffset in [(0, 0), (6, 4)]:
        h[:, offset:offset+4, offset:offset+4] = old[:, oldoffset:oldoffset+4, oldoffset:oldoffset+4]
        v[:, offset:offset+4, offset:offset+4] = vold[:, oldoffset:oldoffset+4, oldoffset:oldoffset+4]
    sl = slice(2*site, 2*site+2)
    for source, target, extra_target in [(old, low, extra), (vold, dlow, dextra)]:
        t = source[:, :4, 4:]
        target[:, :4, 6:10] = t
        target[:, 6:10, :4] = adj(t)
        extra_target[:, :4, 10:12] = t[:, :, sl]
        extra_target[:, 10:12, :4] = adj(t[:, :, sl])
        extra_target[:, 4:6, 6:10] = t[:, sl, :]
        extra_target[:, 6:10, 4:6] = adj(t[:, sl, :])
    pextra = np.diag([0]*4+[1]*2+[0]*4+[1]*2)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    spin = np.kron(np.eye(6), direction[0]*sx+direction[1]*sy)
    bottom = np.diag([0]*6+[1]*6)
    spins = [bottom@spin, (np.eye(12)-bottom)@spin, spin]
    sz = np.kron(np.eye(6), np.diag([1., -1.]))
    return h, v, low, dlow, extra, dextra, pextra, spins, sz, unit.mono.g_factor


def assemble(data, energy, t_low, t_extra, channel='all'):
    h, v, low, dlow, extra, dextra, pextra, _, _, _ = data
    hh = h + energy*pextra + t_low*low + t_extra*extra
    vv = v + t_low*dlow + t_extra*dextra
    if isinstance(channel, dict) or channel != 'all':
        intra, parts, _, _, vi, derivatives = split_channels(
            hh, [np.arange(6), np.arange(6, 12)], 2, vv)
        coefficients = channel if isinstance(channel, dict) else {
            x: 1. for x in (channel.split('_') if channel != 'none' else [])}
        hh = intra + sum((a*parts[x] for x, a in coefficients.items()), np.zeros_like(h))
        vv = vi + sum((a*derivatives[x] for x, a in coefficients.items()), np.zeros_like(v))
    return hh, vv


def selective_extra(data, energy, t_low, t_extra, empty_ratio=0.):
    """Retain the original KM interlayer block; select only new-orbital paths.

    The added valence-empty channel has strength t_extra and its empty-empty
    counterpart has relative amplitude empty_ratio. This is a projected
    effective coupling, generally longer range than the initial orbital hops.
    """
    h, v = assemble(data, energy, 0., t_extra, {'VC': 1., 'CC': empty_ratio})
    return h+t_low*data[2], v+t_low*data[3]


def flatten(h, v, gap=.03):
    """Uniform occupied/empty spectrum, physical readout and Zeeman unchanged."""
    e, u = np.linalg.eigh(h)
    if np.min(e[:, 4]-e[:, 3]) < 1e-9:
        raise ValueError('Cannot flatten a closed gap')
    p = u[:, :, :4]@adj(u[:, :, :4])
    ve = adj(u)@v@u
    dp = np.zeros_like(ve)
    dp[:, :4, 4:] = ve[:, :4, 4:]/(e[:, :4, None]-e[:, None, 4:])
    dp[:, 4:, :4] = adj(dp[:, :4, 4:])
    return gap*(np.eye(h.shape[-1])/2-p), -gap*(u@dp@adj(u))


def response(h, v, spins, sz, weights, step=.002, g=2.):
    """Fixed-four-occupied-band response and band edges for chunk reduction.

    The caller must reject C if the global minimum empty energy is not above
    the global maximum occupied energy. A metallic fixed-filling value is not
    the insulating NMEE and must never be interpreted as an enhancement.
    """
    c = np.zeros(len(spins)); vmax, cmin = -np.inf, np.inf
    for sign in [1, -1]:
        e, u = np.linalg.eigh(h-sign*step*g*MU_B*sz)
        vmax, cmin = max(vmax, e[:, 3].max()), min(cmin, e[:, 4].min())
        occ, empty = u[:, :, :4], u[:, :, 4:]
        vcv = (adj(empty)@v@occ).swapaxes(1, 2)
        den = (e[:, :4, None]-e[:, None, 4:])**2
        for i, s in enumerate(spins):
            pairs = -2*((adj(occ)@s@empty)*vcv).imag/den
            c[i] += sign*(weights@pairs.sum(axis=(1, 2)))/(2*step)
    return c, vmax, cmin
