"""Exact band-sector audit of the gapped T=0 static NMEE; NOT a SW reduction.

The k-dependent middle-four exact eigenspace is separated into two occupied
and two empty states.  The four complementary states must remain in the
response sum unless a controlled dressed-operator reduction is provided.
Orbital embedding and all velocities come from the ORIGINAL 8-band model.
"""
from dataclasses import replace
import numpy as np
from .sliding_tight_binding import SlidingBilayer, SlidingSpec
from nmee_lab.clean_core.core import MU_B

CATEGORY_NAMES = ("low_to_low", "remote_occ_to_low_empty",
                  "low_occ_to_remote_empty", "remote_to_remote")


def band_sector_kernel(H, vx, Sx, Sz, mu, *, g_factor=10.0, step_T=.25):
    """Return [N,4] of d_B chi for four disjoint occ/empty sectors."""
    if step_T <= 0:
        raise ValueError("positive finite Zeeman step required")
    H = np.asarray(H)
    if H.ndim != 3 or H.shape[1:] != (8, 8):
        raise ValueError("expected [N,8,8] Hamiltonian")
    if not np.all((np.linalg.eigvalsh(H)[:, 3] < mu)
                  & (mu < np.linalg.eigvalsh(H)[:, 4])):
        raise ValueError("mu must stay inside the gap")
    field_curvatures = []
    for sign in (1, -1):
        e, U = np.linalg.eigh(H-sign*g_factor*MU_B*step_T*Sz)
        if not np.all((e[:, 3] < mu) & (mu < e[:, 4])):
            raise ValueError("finite Zeeman field closes chemical-potential gap")
        dagger = U.conj().swapaxes(1, 2)
        spin = dagger @ Sx @ U
        velocity = dagger @ vx @ U
        numerator = spin[:, :4, 4:] * np.swapaxes(velocity[:, 4:, :4], 1, 2)
        denom = (e[:, :4, None]-e[:, None, 4:])**2
        pair = -2*np.imag(numerator/denom)
        field_curvatures.append(np.stack((
            pair[:, 2:4, 0:2].sum((1, 2)),
            pair[:, 0:2, 0:2].sum((1, 2)),
            pair[:, 2:4, 2:4].sum((1, 2)),
            pair[:, 0:2, 2:4].sum((1, 2))), axis=-1))
    return (field_curvatures[0]-field_curvatures[1])/(2*step_T)


def valley_sector_integrals(spec: SlidingSpec, *, radius=.5,
                            radial_n=26, angular_n=72, step_T=.25):
    """Independent polar Gauss integration on two K/K' disks (nm^-1).

    Outputs signed contributions, not absolute-integrand capture and not an
    integration bound for the full BZ. Exact spectral bands are used only for
    diagnosis; low_to_low is NOT a complete downfolded response.
    """
    if not radius > 0 or radial_n < 2 or angular_n < 4:
        raise ValueError("invalid disk quadrature")
    model = SlidingBilayer(spec)
    reciprocal = model.model().reciprocal
    centers = np.array([[1/3, 2/3], [2/3, 1/3]]) @ reciprocal
    node, gw = np.polynomial.legendre.leggauss(radial_n)
    r = (node+1)*radius/2
    weights_r = gw*radius/2
    phi = 2*np.pi*(np.arange(angular_n)+.5)/angular_n
    relative = np.stack(np.broadcast_arrays(
        r[:, None]*np.cos(phi)[None, :],
        r[:, None]*np.sin(phi)[None, :]), axis=-1).reshape(-1, 2)
    weights = (np.repeat(weights_r*r, angular_n) *
               (2*np.pi/angular_n)/(2*np.pi)**2)
    val_results, sep, midgap, overlaps = [], [], [], []
    decoupled = SlidingBilayer(replace(spec, tunneling_eV=0.))
    for center in centers:
        k = relative + center
        h = model.matrices(k)
        vx = model.matrices(k, axis=0)
        e, frame = np.linalg.eigh(h)
        e0, frame0 = np.linalg.eigh(decoupled.matrices(k))
        sv = np.linalg.svd(
            frame0[:, :, 2:6].conj().swapaxes(1, 2) @ frame[:, :, 2:6],
            compute_uv=False)
        overlaps.append(float(np.min(sv[:, -1]**2)))
        sep.append(float(min(np.min(e[:, 2]-e[:, 1]),
                             np.min(e[:, 6]-e[:, 5]))))
        midgap.append(float(min(np.min(e[:, 4]-spec.mu_eV),
                                np.min(spec.mu_eV-e[:, 3]))))
        curv = band_sector_kernel(h, vx, model.spins["total_x"],
                                  model.spins["total_z"], spec.mu_eV,
                                  g_factor=model.mono.g_factor, step_T=step_T)
        val_results.append(np.einsum("nc,n->c", curv, weights))
    vals = np.asarray(val_results)
    parts = vals.sum(axis=0)
    return {
        "u_in_shift": list(spec.shift),
        "t_perp_eV": spec.tunneling_eV,
        "radius_inv_nm": radius, "radial_n": radial_n,
        "angular_n": angular_n, "B_step_T": step_T,
        "categories": list(CATEGORY_NAMES),
        "valley_integrals": vals.tolist(), "sum_integrals": parts.tolist(),
        "full": float(parts.sum()), "middle_four_only": float(parts[0]),
        "remote_virtual_transitions": float(parts[1:].sum()),
        "sampled_min_middle4_remote_gap_eV": min(sep),
        "sampled_min_mu_distance_eV": min(midgap),
        "sampled_min_cos2_decoupled_vs_exact_middle4": min(overlaps),
        "protocol": "0K gapped equilibrium static finite-Zeeman response",
    }
