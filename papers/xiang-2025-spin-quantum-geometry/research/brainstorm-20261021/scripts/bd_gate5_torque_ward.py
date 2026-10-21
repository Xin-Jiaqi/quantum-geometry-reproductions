"""B+D Gate 5: layer-odd spin density and interlayer/local torque Ward test.

Purely synthetic 32-band real-space Slater-Koster + atomic SOC + AFM toy.
NOT a reproduction of Xiang et al. (2025), a real material, a new torque
mechanism, or a measured spin transfer.

Run from this scripts/ directory:
  python bd_gate5_torque_ward.py
  python bd_gate5_torque_ward.py --dense
Requires numpy and the neighboring Gate 4 scripts.
Frequency in synthetic energy units; eta=0, below the direct gap.
"""
import json
import sys
import numpy as np
from bd_gate4b_atomic_spin import model, STAG_x, Tspin_x
from bd_gate4_polar_honeycomb import B1, B2, fraction_delta


def response_all(frac=(1/3, 1/3), N=11, omega=.06, staggered=True,
                 m=.34, soc=.30, out_scale=1.):
    """(spin density, interlayer transfer, local torque, sum), minimum direct gap.

    T_inter=i[H_inter,S]; T_local=i[H-H_inter,S].
    Retarded complete occupied-unoccupied length-gauge response at eta=0.
    """
    S = STAG_x if staggered else Tspin_x
    delta = fraction_delta(frac)
    vals = np.zeros(4, complex)
    direct_gap = np.inf
    for i in range(N):
        for j in range(N):
            k = (i * B1 + j * B2) / N
            H, V = model(k, delta, m=m, soc=soc, out_scale=out_scale)
            H0, _ = model(k, delta, m=m, soc=soc, out_scale=0.)
            Hi = H - H0
            transfer = 1j * (Hi @ S - S @ Hi)
            local = 1j * (H0 @ S - S @ H0)
            total = 1j * (H @ S - S @ H)
            assert np.linalg.norm(transfer + local - total) < 2e-13
            E, U = np.linalg.eigh(H)
            direct_gap = min(direct_gap, E[8] - E[7])
            V_band = U.conj().T @ V @ U
            gap = E[8:][None, :] - E[:8, None]
            for n, O in enumerate((S, transfer, local, total)):
                O_band = U.conj().T @ O @ U
                q = O_band[:8, 8:] * V_band[8:, :8].T
                r = -1j * q / gap
                vals[n] += np.sum(
                    r / (omega - gap) - r.conj() / (omega + gap)
                ) / N**2
    return vals, direct_gap


def check(N, omega=.06):
    ab, gap = response_all(N=N, omega=omega)
    ba, gap_ba = response_all(frac=(-1/3, -1/3), N=N, omega=omega)
    assert gap > omega and gap_ba > omega
    assert np.max(np.abs(ab + ba)) < 8e-13  # physical m_z T partner
    assert abs(ab[1] + ab[2] - ab[3]) < 1e-12
    assert abs(ab[3] + 1j * omega * ab[0]) < 1e-12
    assert abs(ab[1]) > 5e-5 and abs(ab[3]/ab[1]) < .02

    total, _ = response_all(N=N, omega=omega, staggered=False)
    # The interlayer hopping of the full 32-band model is spin independent.
    assert abs(total[1]) < 1e-13

    for kw in (dict(out_scale=0.), dict(soc=0.), dict(m=0.)):
        off, _ = response_all(N=N, omega=omega, **kw)
        assert max(abs(v) for v in off) < 1e-11

    return {
        "mesh": [N, N], "omega": omega, "eta": 0.,
        "minimum_direct_gap": float(gap),
        "chi_S_staggered_x": [ab[0].real, ab[0].imag],
        "chi_T_inter": [ab[1].real, ab[1].imag],
        "chi_T_local": [ab[2].real, ab[2].imag],
        "chi_T_total": [ab[3].real, ab[3].imag],
        "ward_residual": float(abs(ab[3] + 1j*omega*ab[0])),
        "remaining_ratio": float(abs(ab[3]/ab[1])),
        "switch_parity_residual": float(np.max(np.abs(ab+ba))),
        "total_spin_interlayer_torque_abs": float(abs(total[1])),
    }


if __name__ == "__main__":
    dense = "--dense" in sys.argv
    for n in ([11, 19, 23] if dense else [7, 11]):
        print(json.dumps(check(n), indent=2))
    print("MODEL PASS; NOVELTY UNRESOLVED; NO REAL MATERIAL CLAIM")
