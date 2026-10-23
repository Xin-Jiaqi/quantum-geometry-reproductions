"""Stage8n Kubo/projector static-response core (local, testable).

F(B) = sum_{n occ, m unocc} 2 Im[ <m|S_i/2|n> <n|v|m> ] / (E_n - E_m)^2
with <m|v|n> = (E_n - E_m) <m|d u_n/dk> (off-diagonal part of dH/dk in the
band basis). All operators are in one common basis; the off-diagonal
Im[S_mn v_mn] is gauge invariant under unitary column-phase rotations of
that basis (used by the gauge-covariance regression test below).
"""
import numpy as np


def kubo_gamma_k(E, S2, dM, occ):
    """Gamma_k from energies, spin source S2 (=SPN/2 in band basis) and
    band-basis overlap gradient dM[m,n] = <u_m|d u_n/dk> (central diff of
    MMN overlaps along the k axis).

    From the resolvent expansion:
      F = +2 sum_{n occ, m unocc} Im[ <n|S_i/2|m> <m|d u_n/dk> ] / (E_n - E_m)
    (the velocity matrix element <m|dH/dk|n> = (E_n-E_m)<m|d u_n/dk> cancels
    one power of the denominator; the <n|S|m> index order is the conjugate
    of <m|S|n>, as required by the a* b pairing of the P[.,.] trace).

    Returns (Gamma_k, min occ-unocc gap).
    """
    nb = len(E)
    gam = 0.0
    gmin = float("inf")
    for n in range(occ):
        for m in range(occ, nb):
            d = E[n] - E[m]
            gmin = min(gmin, abs(d))
            if abs(d) > 1e-6:
                snm = S2[n, m]             # <n|S_i/2|m>
                dumn = dM[m, n]            # <m|d u_n/dk>
                gam += 2.0 * np.imag(snm * dumn) / d
    return float(gam), float(gmin)


def rotate_to_eigenbasis(dM, S2, U):
    """Rotate both operators by the SAME unitary U (Zeeman eigenbasis):
    dM_b = U^+ dM U, S2_b = U^+ S2 U. Gauge-consistent construction used
    by the real-data diagnostic (diag3)."""
    return (U.T.conj() @ dM @ U,
            U.T.conj() @ S2 @ U)
