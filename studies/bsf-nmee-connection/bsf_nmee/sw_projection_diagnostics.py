"""Diagnostic ONLY: test a fixed four-band valley projector as a SW starting point.

P0 consists of the middle four eigenvectors of the decoupled bilayer at K/K'.
This is not a downfolded NMEE formula, and ||QHP||/Delta is not a rigorous
Schrieffer-Wolff remainder bound. No frozen paper reproduction is modified.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.models import SY
from .sliding_tight_binding import SlidingBilayer, SlidingSpec, M_GLOBAL

TR_UNITARY = np.kron(np.eye(4), 1j * SY)


def sw_gate(*, u_values=(.1, .4), tunneling_values=(.005, .02),
            q_x_values=(0., .1, .25, .5, .75)):
    """Return 40 rows for the default parameters, two valleys and a radial line.

    k is in nm^-1; q=(q_x,0). A full disk or BZ bound is NOT established.
    Four-band subspace is fixed at t_perp=0, q=0 for each valley.
    """
    base = SlidingBilayer(replace(SlidingSpec(), shift=(0., 0.),
                                   tunneling_eV=0.))
    g = base.model().reciprocal
    centers = np.array([[1/3, 2/3], [2/3, 1/3]]) @ g
    output = []
    for iv, center in enumerate(centers):
        _, u0 = np.linalg.eigh(base.matrices(center[None, :])[0])
        p = u0[:, 2:6]
        q = np.concatenate([u0[:, :2], u0[:, 6:]], axis=1)
        for us in u_values:
            for t in tunneling_values:
                shift = tuple(us * np.array([.23, .41]))
                model = SlidingBilayer(replace(SlidingSpec(),
                                               shift=shift, tunneling_eV=t))
                mirror_model = SlidingBilayer(replace(SlidingSpec(),
                                                      shift=tuple(-np.array(shift)),
                                                      tunneling_eV=t))
                for radius in q_x_values:
                    k = center + np.array([radius, 0.])
                    h = model.matrices(k[None, :])[0]
                    a = p.conj().T @ h @ p
                    d = q.conj().T @ h @ q
                    b = q.conj().T @ h @ p
                    ea = np.linalg.eigvalsh(a)
                    ed = np.linalg.eigvalsh(d)
                    delta = float(np.min(np.abs(ea[:, None] - ed[None, :])))
                    norm = float(np.linalg.norm(b, ord=2))
                    e, eigvec = np.linalg.eigh(h)
                    sv = np.linalg.svd(p.conj().T @ eigvec[:, 2:6],
                                       compute_uv=False)
                    hm = mirror_model.matrices(k[None, :])[0]
                    hminus = model.matrices((-k)[None, :])[0]
                    mirror_err = float(np.max(np.abs(h - M_GLOBAL @ hm @ M_GLOBAL.conj().T)))
                    tr_err = float(np.max(np.abs(h - TR_UNITARY @ hminus.conj() @ TR_UNITARY.conj().T)))
                    output.append({
                        "valley": iv, "u": float(us),
                        "t_perp_meV": float(t * 1000),
                        "q_x_inv_nm": float(radius),
                        "p_q_coupling_meV": norm * 1000,
                        "block_min_spectral_distance_meV": delta * 1000,
                        "sw_diagnostic_ratio": norm / delta if delta else float("inf"),
                        "full_middle_gap_meV": float((e[4] - e[3]) * 1000),
                        "mu_distance_meV": float(min(e[4]-model.spec.mu_eV,
                                                    model.spec.mu_eV-e[3]) * 1000),
                        "smallest_principal_cos_squared": float(min(sv**2)),
                        "mirror_error_eV": mirror_err, "TR_error_eV": tr_err,
                    })
    return output
