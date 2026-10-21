"""v64: strictly unperturbed two-drive-contact limit of virtual QWZ voltage probes.

The voltage is the lambda -> 0+ limit of floating TOP/BOTTOM probes,
calculated from a retarded *two-lead* Green function with no probe
selfenergies. This mathematical limit is different from measuring at
lambda=0, where the two floating voltages are undetermined. Do not regard
finite-L, finite-eta results as a thermodynamic Hall resistivity, or as an
orbital mixed E/B response.
"""
import numpy as np
from .qwz_finite_marker_v60 import qwz_open
from .qwz_spatial_v62 import qwz_width_leads
from .matrix_leads_v60 import selfenergy_matrix


def virtual_probe_limit(L=7, mass=-1., drive_width=3, probe_width=3,
                        mu=0., eta=1e-4, offsets=(0, 0, 0, 0),
                        min_probe_transmission=1e-13, min_drive_current=1e-7):
    """Compute lambda->0+ floating voltages from independent two-lead G0.

    Vp = [T_pL VL + T_pR VR]/[T_pL + T_pR]; every T_p* is evaluated
    with G0 and a *unit amplitude* probe gamma, before multiplying by lambda².
    This ratio is well-defined only for nonzero probe spectral coupling;
    the completely detached probe's voltage is never claimed measurable.
    Rxy = (Vtop - Vbottom) / two-lead I_L, with V_L=+1/2,V_R=-1/2.
    """
    if not isinstance(L, int) or L < 3:
        raise ValueError('integer L >= 3')
    if not isinstance(drive_width, int) or not 2 <= drive_width <= L:
        raise ValueError('drive width 2..L')
    if not isinstance(probe_width, int) or not 2 <= probe_width <= L:
        raise ValueError('probe width 2..L')
    if len(offsets) != 4 or any(not isinstance(x, int) for x in offsets):
        raise ValueError('four integer offsets')
    if not np.isfinite([mass, mu, eta, min_probe_transmission, min_drive_current]).all():
        raise ValueError('finite numeric inputs')
    if eta <= 0 or min_probe_transmission <= 0 or min_drive_current <= 0:
        raise ValueError('positive eta and measurement thresholds')
    H, _ = qwz_open(L, L, mass)
    drive_leads = qwz_width_leads(L, L, mass, drive_width, eta,
                                 (offsets[0], 0, offsets[2], 0))
    probe_leads = qwz_width_leads(L, L, mass, probe_width, eta,
                                 (0, offsets[1], 0, offsets[3]))
    left, right = drive_leads[0], drive_leads[2]
    top, bottom = probe_leads[1], probe_leads[3]
    S_L, _, G_L, _ = selfenergy_matrix(mu, left)
    S_R, _, G_R, _ = selfenergy_matrix(mu, right)
    _, _, G_T, _ = selfenergy_matrix(mu, top)
    _, _, G_B, _ = selfenergy_matrix(mu, bottom)
    Gr = np.linalg.inv(mu*np.eye(len(H)) - H - S_L - S_R)
    Ga = Gr.conj().T
    a_left = Gr@G_L@Ga
    a_right = Gr@G_R@Ga
    # Two-terminal Caroli is completely independent of floating probes.
    T_lr = float(np.trace(G_R@a_left).real)
    T_rl = float(np.trace(G_L@a_right).real)
    weights = np.array([[np.trace(g@a_left).real, np.trace(g@a_right).real]
                        for g in (G_T,G_B)], dtype=float)
    if np.min(weights) < -1e-10:
        raise RuntimeError('negative virtual transmission')
    weights = np.maximum(0., weights)
    totals = weights.sum(axis=1)
    measurable = np.all(totals > min_probe_transmission) and abs(T_lr) >= min_drive_current
    voltages = (weights[:,0]-weights[:,1])/(2*totals) if np.all(totals > min_probe_transmission) else None
    rxy = float((voltages[0]-voltages[1])/T_lr) if measurable else None
    return {
        'L':L,'mass':float(mass),'drive_width':drive_width,'probe_width':probe_width,
        'mu':float(mu),'eta':float(eta),'offsets':tuple(offsets),
        'two_terminal_T_lr':T_lr,'two_terminal_T_rl':T_rl,
        'probe_spectral_weights':weights, 'probe_weight_totals':totals,
        'floating_voltages_lambda_to_0':voltages,
        'Rxy_lambda_to_0_h_over_q2':rxy,
        'limit_defined':bool(measurable),
        'spectral_identity_error':float(np.linalg.norm(1j*(Gr-Ga)-a_left-a_right)),
        'status':'VIRTUAL_LAMBDA_TO_ZERO_PLUS_FINITE_QWZ_DEVICE_NOT_ZERO_COUPLING_MEASUREMENT_NOT_BULK_EB'
    }
