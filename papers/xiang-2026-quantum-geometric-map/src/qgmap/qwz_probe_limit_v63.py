"""v63: non-invasive voltage-probe limit for the same microscopic QWZ device.

Only the *device-to-probe coupling amplitude* of top/bottom voltage probes is
scaled by lambda. Their lead H00/H01 and eta are unchanged; Sigma and Gamma
therefore scale as lambda**2. Left/right drive electrodes remain fixed.

This is a coherent finite-device probe-limit study, not a thermodynamic
Kubo--NEGF equivalence proof and not the original orbital mixed E/B response.
"""
import numpy as np
from .qwz_finite_marker_v60 import qwz_open
from .qwz_spatial_v62 import qwz_width_leads
from .matrix_leads_v60 import MatrixLead, selfenergy_matrix
from .lead_transport_v57 import four_terminal


def probe_strength_response(L=7, mass=-1., width=3, mu=0., eta=1e-4,
                            probe_scale=1e-3, offsets=(0,0,0,0), probe_width=None):
    """Return four-probe transport, local bond currents, and lead broadening.

    lambda=probe_scale multiplies C for top and bottom only. Since
    Sigma=C g C†, Gamma scales lambda² exactly up to numerical precision.
    Drive-contact width and probe-contact width are independent controls.
    lambda must be strictly positive: lambda=0 is not a voltage measurement
    and makes the floating-probe voltage equations singular.
    """
    if not isinstance(L, int) or L < 3:
        raise ValueError("L must be integer >=3")
    if probe_width is None:
        probe_width = width
    if not isinstance(width, int) or not 2 <= width <= L:
        raise ValueError("drive width must be integer 2..L")
    if not isinstance(probe_width, int) or not 2 <= probe_width <= L:
        raise ValueError("probe_width must be integer 2..L")
    if not np.isfinite([mass, mu, eta, probe_scale]).all() or eta <= 0 or probe_scale <= 0:
        raise ValueError("finite mass/mu, eta>0 and probe_scale>0 required")
    if len(offsets) != 4 or any(not isinstance(x, int) for x in offsets):
        raise ValueError("four integer side offsets")

    H, xy = qwz_open(L, L, mass)
    # Drive-contact width and probe-contact width are independent controls.
    drive_offsets = (offsets[0], 0, offsets[2], 0)
    probe_offsets = (0, offsets[1], 0, offsets[3])
    drive_leads = qwz_width_leads(L, L, mass, width, eta, drive_offsets)
    probe_leads = qwz_width_leads(L, L, mass, probe_width, eta, probe_offsets)
    left, right = drive_leads[0], drive_leads[2]
    top, bottom = probe_leads[1], probe_leads[3]
    leads = (
        left,
        MatrixLead(top.h00, top.h01, probe_scale*top.coupling, top.eta),
        right,
        MatrixLead(bottom.h00, bottom.h01, probe_scale*bottom.coupling, bottom.eta),
    )

    parts = [selfenergy_matrix(mu, lead) for lead in leads]
    Sigma = sum((p[0] for p in parts), np.zeros_like(H))
    Gammas = [p[2] for p in parts]
    Gr = np.linalg.inv(mu*np.eye(len(H)) - H - Sigma)
    Ga = Gr.conj().T
    injected_spectral = [Gr @ G @ Ga for G in Gammas]

    T = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if a != b:
                T[a, b] = np.trace(Gammas[a] @ injected_spectral[b]).real
    obs = four_terminal(T, np.zeros_like(T))
    rho_bias = sum((v*a for v, a in zip(obs["voltages"], injected_spectral)),
                   np.zeros_like(H))

    A = 1j*(Gr-Ga)
    spectral_error = float(np.linalg.norm(A-sum(injected_spectral)))
    ldos = np.real(np.diag(A)).reshape(L, L, 2).sum(axis=2)/(2*np.pi)
    Jx = np.zeros((L, L-1))
    Jy = np.zeros((L-1, L))

    def bond_current(i, j):
        si, sj = slice(2*i, 2*i+2), slice(2*j, 2*j+2)
        return float(-2*np.imag(np.trace(H[si, sj] @ rho_bias[sj, si])))

    for y in range(L):
        for x in range(L):
            i = y*L+x
            if x < L-1:
                Jx[y, x] = bond_current(i, i+1)
            if y < L-1:
                Jy[y, x] = bond_current(i, i+L)

    divergence = np.zeros((L, L))
    for y in range(L):
        for x in range(L):
            divergence[y, x] = ((Jx[y, x] if x < L-1 else 0)
                                -(Jx[y, x-1] if x else 0)
                                +(Jy[y, x] if y < L-1 else 0)
                                -(Jy[y-1, x] if y else 0))

    IL = float(obs["longitudinal_conductance_q2_over_h"])
    reportable = abs(IL) >= 1e-7
    return {
        "L": L, "mass": float(mass), "width": width, "probe_width": probe_width,
        "mu": float(mu), "eta": float(eta), "probe_scale": float(probe_scale),
        "transmissions": T, "voltages": obs["voltages"],
        "drive_current_q2_over_h": IL,
        "Rxy_h_over_q2": float(obs["hall_ratio"]/IL) if reportable else None,
        "Rxy_reportable": bool(reportable),
        "Jx": Jx, "Jy": Jy, "bond_divergence": divergence,
        "spectral_ldos": ldos, "spectral_identity_error": spectral_error,
        "x_cut_currents": Jx.sum(axis=0), "y_cut_currents": Jy.sum(axis=1),
        "floating_current_error": obs["max_floating_probe_current"],
        "total_lead_current_error": obs["max_current_conservation_residual"],
        "coupling_norms": np.array([np.linalg.norm(l.coupling) for l in leads]),
        "gamma_norms": np.array([np.linalg.norm(g) for g in Gammas]),
        "selfenergy_norms": np.array([np.linalg.norm(p[0]) for p in parts]),
        "probe_row_transmission": float(T[1, :].sum()+T[3, :].sum()),
        "max_interior_divergence": float(np.max(np.abs(divergence[1:-1, 1:-1]))),
        "status": "CONDITIONAL_FINITE_DEVICE_NONINVASIVE_VOLTAGE_PROBE_LIMIT; NOT_ORBITAL_EB",
    }
