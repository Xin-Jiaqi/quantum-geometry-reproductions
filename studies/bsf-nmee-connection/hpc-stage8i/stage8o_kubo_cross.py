"""Stage8o corrected Kubo/projector static-response cross-check.

Two independent estimators of
  F(B) = i Tr P[ dP/dlam , dP/dk ]   (lambda=source spin axis, k=axis)
in the frozen-KS frame H(k,lam,B) = diag E + lam*S_i/2 - muB*B*S_z:

(1) RESOLVENT (analytic Kubo):
      F = 2 sum_{n occ, m unocc} Im[ S2[n,m] dM[m,n] ] / (E_n - E_m)
  with S2 = <n|S_i/2|m> in the Zeeman eigenbasis at B and dM the
  gauge-covariant overlap gradient from the MMN links,
      dM = (L(k,k+dk) - L(k,k-dk))/(2 dk),  L(k,k') = U(B,k)^+ M(k,k') U(B,k')
  (the transported links include the U^+ d_k U term Stage8n dropped).

(2) PROJECTOR finite differences (independent of perturbation theory):
      dP/dlam ~ [P(lam+dl) - P(lam-dl)]/(2 dl)  in the u(B,k) basis
      dP/dk   ~ dM  (band-basis matrix of dP/dk is exactly dM for m!=n)
      F = i Tr[ P (dP/dlam)(dP/dk) - P (dP/dk)(dP/dlam) ]
  P is the diagonal occ-block in the u(B,k) basis.

Both must agree on a synthetic nonzero, well-gapped model before real use.
"""
import numpy as np


def dM_covariant(U, M_plus, M_minus, dk):
    """Gauge-covariant band-basis overlap gradient incl. U^+ d_k U:
    dM[m,n] = <u_m(B,k)|d u_n(B,k)/dk> from two transported links."""
    Lp = U.T.conj() @ M_plus @ U
    Lm = U.T.conj() @ M_minus @ U
    return (Lp - Lm) / (2.0 * dk)


def resolvent_gamma_k(E, S2_zb, dM, occ):
    """(1) Analytic Kubo: 2 Im[S2[n,m] dM[m,n]]/(E_n-E_m), n occ, m unocc."""
    nb = len(E)
    gam = 0.0
    gmin = float("inf")
    for n in range(occ):
        for m in range(occ, nb):
            d = E[n] - E[m]
            gmin = min(gmin, abs(d))
            if abs(d) > 1e-6:
                gam += 2.0 * np.imag(S2_zb[n, m] * dM[m, n]) / d
    return float(gam), float(gmin)


def projector_gamma_k(E, Sz_vasp, Sx_vasp, M_plus, M_minus, occ, dk, dl, B,
                      mu_b=5.7884e-5, g=2.0):
    """(2) Direct iTr P[dP/dlam, dP/dk] by finite differences.

    All matrices are (nb, nb) in the VASP band basis; the Zeeman basis at B
    is U_B = eigh(H(B)). The projector-route derivative dP/dk is taken to be
    the band-basis matrix dM (gauge-covariant transported links), which is
    exact for the off-diagonal (occ-unocc) part that the trace uses.
    """
    H_B = np.diag(E).astype(complex) - g * mu_b * B * Sz_vasp
    e_B, U_B = np.linalg.eigh(H_B)
    P_occ = np.zeros_like(U_B)
    P_occ[:occ, :occ] = 1.0
    # dP/dlam in the u(B) basis: projectors at lam +/- dl rebuilt by
    # re-diagonalizing H(lam +/- dl, B) and transported to u(B).
    def P_of_lam(dlam):
        H = H_B + dlam * Sx_vasp
        _, V = np.linalg.eigh(H)
        Pv = V[:, :occ] @ V[:, :occ].T.conj()
        return U_B.T.conj() @ Pv @ U_B
    dPdl = (P_of_lam(dl) - P_of_lam(-dl)) / (2.0 * dl)
    dM = dM_covariant(U_B, M_plus, M_minus, dk)
    F = 1j * np.trace(P_occ @ dPdl @ dM - P_occ @ dM @ dPdl)
    return float(F.real), float(np.linalg.norm(dPdl, "fro"))


def zeeman_energies(E, Sz_vasp, B, mu_b=5.7884e-5, g=2.0):
    H = np.diag(E).astype(complex) - g * mu_b * B * Sz_vasp
    return np.linalg.eigh(H)[0]
