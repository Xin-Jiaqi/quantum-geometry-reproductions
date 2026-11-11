"""Stage 7c: exact channel-sign / spectral-sector NMEE audit at K/K'.

Physical scope: eight-band, 0-K, insulating equilibrium STATIC finite-Zeeman
spin-electric response. Not driven DC, Wannier DFT, or a ferroelectric state.
The declared same/cross orbital-channel signs are basis-dependent model probes.
"""
from dataclasses import replace
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from .spectral_sector import band_sector_kernel, CATEGORY_NAMES


def disk_channel_sectors(spec: SlidingSpec, *, same_scale=1., cross_scale=1.,
                         radius=.75, radial_n=26, angular_n=72, step_T=.25):
    """Exact four occupied-to-empty band sectors on two nonoverlapping disks."""
    if radius <= 0 or radial_n < 2 or angular_n < 4 or step_T <= 0:
        raise ValueError("invalid quadrature radius, counts or field step")
    obj = SlidingBilayer(spec)
    obj.hops = [(ia, ib, d, amp*(same_scale if ia == ib else cross_scale))
                for ia, ib, d, amp in obj.hops]
    reciprocal = obj.model().reciprocal
    centers = np.array([[1./3, 2./3], [2./3, 1./3]]) @ reciprocal
    nodes, w_rad = np.polynomial.legendre.leggauss(radial_n)
    rr = (nodes+1)*radius/2
    dr = w_rad*radius/2
    phi = 2*np.pi*(np.arange(angular_n)+.5)/angular_n
    displacement = (rr[:, None, None] * np.stack((
        np.cos(phi), np.sin(phi)), axis=-1)[None, :, :]).reshape(-1, 2)
    weights = np.repeat(rr*dr, angular_n)/(2*np.pi*angular_n)
    spin_x = obj.spins["total_x"]
    spin_z = obj.spins["total_z"]
    sectors = []
    mu_gaps = []
    outer_gaps = []
    for center in centers:
        k = displacement+center
        H = obj.matrices(k)
        Vx = obj.matrices(k, axis=0)
        eig = np.linalg.eigvalsh(H)
        mu_gap = float(min(np.min(eig[:,4]-spec.mu_eV),
                           np.min(spec.mu_eV-eig[:,3])))
        if mu_gap <= 0:
            raise ValueError("not insulating at chosen chemical potential")
        mu_gaps.append(mu_gap)
        outer_gaps.append(float(min(np.min(eig[:,2]-eig[:,1]),
                                     np.min(eig[:,6]-eig[:,5]))))
        kernel = band_sector_kernel(
            H, Vx, spin_x, spin_z, spec.mu_eV,
            g_factor=obj.mono.g_factor, step_T=step_T)
        sectors.append(np.einsum("nc,n->c", kernel, weights))
    sectors = np.asarray(sectors)
    parts = sectors.sum(axis=0)
    return {
        "valley_sectors": sectors.tolist(),
        "sectors": parts.tolist(), "total": float(parts.sum()),
        "sampled_min_mu_distance_eV": min(mu_gaps),
        "sampled_min_outer_group_gap_eV": min(outer_gaps),
        "radius_inv_nm": radius, "radial_n": radial_n,
        "angular_n": angular_n, "step_T": step_T,
        "same_scale": same_scale, "cross_scale": cross_scale,
        "categories": list(CATEGORY_NAMES),
    }


def channel_mechanism(spec: SlidingSpec, **quad):
    """A, B, and odd-mixed/even coefficients normalized by nominal t_perp².

    Mixed-odd = (Gamma[t,t]-Gamma[t,-t])/(2*t²), not *exactly* the
    weak-coupling C unless t -> 0. The even part similarly contains t⁴ terms.
    """
    t = float(spec.tunneling_eV)
    if not np.isfinite(t) or t <= 0:
        raise ValueError("positive reference tunneling required")
    choices = {"full": (1.,1.), "cross_reverse": (1.,-1.),
               "same_only": (1.,0.), "cross_only": (0.,1.),
               "zero": (0.,0.)}
    cases = {name: disk_channel_sectors(
        spec,same_scale=s,cross_scale=c,**quad)
        for name,(s,c) in choices.items()}
    parts = {name: np.asarray(val["sectors"]) for name,val in cases.items()}
    norm = t*t
    odd = (parts["full"]-parts["cross_reverse"])/(2*norm)
    even = (parts["full"]+parts["cross_reverse"])/(2*norm)
    a = (parts["same_only"]-parts["zero"])/norm
    b = (parts["cross_only"]-parts["zero"])/norm
    return {
        "u_shift_fractional": list(spec.shift), "t_eV": t,
        "mixed_odd_sectors": odd.tolist(), "even_sectors": even.tolist(),
        "same_only_sectors": a.tolist(), "cross_only_sectors": b.tolist(),
        "mixed_total": float(odd.sum()), "even_total": float(even.sum()),
        "same_only_total": float(a.sum()), "cross_only_total": float(b.sum()),
        "full": cases["full"]["total"],
        "full_low_low_over_t2": float(parts["full"][0]/norm),
        "full_remote_over_t2": float(parts["full"][1:].sum()/norm),
        "min_mu_distance_eV": min(v["sampled_min_mu_distance_eV"]
                                 for v in cases.values()),
        "by_case": cases,
        "scope": "two valley DISKS; gapped 0-K equilibrium static, not DC",
    }
