#!/usr/bin/env python3
"""Stage 27 — source-gauge projector geometry vs cell-gauge COVARIANT geometry.

Independent audit requested by the Stage 26 review (item B). All conventions
are strictly unified across every comparison:

  * Hamiltonian:   analytic 2-band spinless Haldane H_src (M=0.4, phi=pi/2,
                   t1=1.0, t2=0.2; NN=(3,2), NNN=(3,2)) and its gauge rotation
                   H_cell(k) = U^dagger(k) H_src(k) U(k), U = diag(1, e^{-ik.d1}),
                   d1=(1,0) A. Cross-checked against the synthetic package HR
                   kernel (G0).
  * Wannier centers: tau_a = i (d_a U) U^dagger (Cartesian k units), verified
                   against the package wannier90-wcenter.dat differences (G0b).
  * Covariant derivative: D_a X = d_a X + i[X, tau_a]; the physically matched
                   object in the cell gauge is D_a P_cell, NOT d_a P_cell.
  * k mesh / broadening / band regulator / tensor conventions: identical to the
                   Stage 26 oracle (5x5x1 (+9,13), eta_sp=0.05, eta_band scan,
                   components xxx/xyy/yxy with (bc) symmetrization, f_nm delta).

Relations established:
  G1  D_a H_cell         == U^dagger (d_a H_src) U
  G2  D_a P_cell         == U^dagger (d_a P_src) U
  G3  sigma(P_src, d)    == sigma(P_cell, D)            (gauge covariance)
  G4  sigma(P_cell, d)   != sigma(P_cell, D)            (plain cell-d is wrong)
  G5  sigma(P_cell, D)   == -pi * R_dressed             (NLO production raw IS the
                                                          covariant response)
  G6  sigma(P_cell, d)   == -pi * R_cell_plain          (Stage 26 plain pair,
                                                          gauge-dependent, kept)
  G7  R_dressed          == R_src                       (production raw == source
                                                          gauge raw)
  G8  sigma(P_src, d)    == -pi * R_src
  G5b eta_band and k-grid stability of G5.
Pairwise table of ALL (sigma, -pi*raw) pairings is reported so mismatched
pairings are documented, not hidden.

The -pi/2 relation claimed in Stage 25 is NOT assumed anywhere; the measured
constant sigma/(-raw) is reported per pairing.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_HERE = Path(__file__).parent
sys.path.insert(0, str(_HERE))

from stage26_oracle import (  # noqa: E402
    COMPONENTS,
    ETA_SP,
    OMEGAS,
    read_centers,
    read_hr,
    reduced_shift,
    fourier_for_sigma,
    grid_kpoints,
    model_raw,
)
from stage26_second_derivative_validation import (  # noqa: E402
    B1,
    B2,
    analytic_source_hamiltonian,
)

SQ3 = np.sqrt(3.0)
B_MAT = np.vstack([B1[:2], B2[:2]])  # reduced -> Cartesian (A^-1), 2D
DELTA1 = np.array([1.0, 0.0])

MODEL_DIR = _HERE.parent / "models" / "stage24_haldane_synthetic" / "L1"

ETA_BAND_PRIMARY = 0.005
H_D1, H_D2 = 1e-4, 3e-3
ETA_BAND_GRID = (0.05, 0.02, 0.01, 0.005, 0.002)
KGRID = (5, 9, 13)
STRONG_RAW = 0.02  # strong-point mask threshold on |raw|
NW = 2


def kcart_of(kpts):
    return np.column_stack([kpts[:, :2] @ B_MAT, np.zeros(len(kpts))])


def h_src_at(kpts):
    out = np.zeros((len(kpts), NW, NW), dtype=np.complex128)
    for i, k in enumerate(kcart_of(kpts)):
        out[i] = analytic_source_hamiltonian(k)[0]
    return out


def d1_d2_src_at(kpts):
    """Analytic source derivatives in the (nk, nw, nw, 3[,3]) convention
    (derivative index LAST, matching NLO/oracle band-matrix layout)."""
    nk = len(kpts)
    dH = np.zeros((nk, NW, NW, 3), dtype=np.complex128)
    d2H = np.zeros((nk, NW, NW, 3, 3), dtype=np.complex128)
    for i, k in enumerate(kcart_of(kpts)):
        _, d1, d2 = analytic_source_hamiltonian(k)
        dH[i] = d1.transpose(1, 2, 0)
        d2H[i] = d2.transpose(2, 3, 0, 1)
    return dH, d2H


def u_at(kpts):
    """U = diag(1, e^{-ik.d1}) in the Cartesian-k convention."""
    nk = len(kpts)
    U = np.zeros((nk, NW, NW), dtype=np.complex128)
    for i, k in enumerate(kcart_of(kpts)):
        U[i] = np.diag([1.0 + 0j, np.exp(-1j * float(k[:2] @ DELTA1))])
    return U


def h_cell_at(kpts):
    """H_cell(k) = U^dagger H_src U on a k grid."""
    U = u_at(kpts)
    Ud = U.conj().transpose(0, 2, 1)
    return np.einsum("kij,kjm,kmo->kio", Ud, h_src_at(kpts), U, optimize=True)


def tau_from_U(kpts, h=H_D1):
    """tau_a = i (d_a U) U^dagger in Cartesian units, via central FD."""
    U0 = u_at(kpts)
    Ud0 = U0.conj().transpose(0, 2, 1)
    tau = np.zeros((3, NW, NW), dtype=np.complex128)
    for a in range(3):
        shift = reduced_shift(h * np.eye(3)[a])
        dU = (u_at(kpts + shift) - u_at(kpts - shift)) / (2 * h)
        tau[a] = np.einsum("kij,kjm->im", dU, Ud0) / len(kpts)
    tau = 1j * tau
    # tau must be the real diagonal center matrix; report the deviation
    return tau


def projectors_from_H(H_at, ks):
    """Eigenstate projectors P[:, :, :, n] on a k grid."""
    nk = len(ks)
    _, wf = np.linalg.eigh(H_at(ks))[:2]
    P = np.zeros((nk, NW, NW, NW), dtype=np.complex128)
    for n in range(NW):
        u = wf[:, :, n]
        P[:, :, :, n] = u[:, :, None] * u[:, None, :].conj()
    return P


def plain_geometry(H_of, kpts, h1=H_D1, h2=H_D2):
    """P, d_a P, d_a d_b P with ORDINARY derivatives (flat-connection gauge)."""
    nk = len(kpts)
    P = projectors_from_H(H_of, kpts)
    dP = np.zeros((nk, NW, NW, NW, 3), dtype=np.complex128)
    for a in range(3):
        sa = reduced_shift(h1 * np.eye(3)[a])
        dP[..., a] = (projectors_from_H(H_of, kpts + sa)
                      - projectors_from_H(H_of, kpts - sa)) / (2 * h1)
    d2P = np.zeros((nk, NW, NW, NW, 3, 3), dtype=np.complex128)
    for a in range(3):
        sa = reduced_shift(h2 * np.eye(3)[a])
        for b in range(3):
            sb = reduced_shift(h2 * np.eye(3)[b])
            d2P[..., a, b] = (
                projectors_from_H(H_of, kpts + sa + sb)
                - projectors_from_H(H_of, kpts + sa - sb)
                - projectors_from_H(H_of, kpts - sa + sb)
                + projectors_from_H(H_of, kpts - sa - sb)
            ) / (4 * h2 * h2)
    return P, dP, d2P


def covariant_geometry(H_of, kpts, tau, h1=H_D1, h2=H_D2):
    """P, D_a P, D_a D_b P with D_a X = d_a X + i[X, tau_a]."""
    nk = len(kpts)

    def D_bP_at(ks, b, h):
        sb = reduced_shift(h * np.eye(3)[b])
        Pk = projectors_from_H(H_of, ks)
        dP = (projectors_from_H(H_of, ks + sb)
              - projectors_from_H(H_of, ks - sb)) / (2 * h)
        return dP + 1j * (np.einsum("kijb,jm->kimb", Pk, tau[b], optimize=True)
                          - np.einsum("ij,kjmb->kimb", tau[b], Pk, optimize=True))

    P = projectors_from_H(H_of, kpts)
    dP = np.zeros((nk, NW, NW, NW, 3), dtype=np.complex128)
    for b in range(3):
        dP[..., b] = D_bP_at(kpts, b, h2)
    d2P = np.zeros((nk, NW, NW, NW, 3, 3), dtype=np.complex128)
    for a in range(3):
        sa = reduced_shift(h1 * np.eye(3)[a])
        for b in range(3):
            # D_a(D_b P) = d_a(D_b P) + i[D_b P, tau_a]
            outer = (D_bP_at(kpts + sa, b, h2) - D_bP_at(kpts - sa, b, h2)) / (2 * h1)
            d2P[..., a, b] = outer + 1j * (
                np.einsum("kijb,jm->kimb", dP[..., b], tau[a], optimize=True)
                - np.einsum("ij,kjmb->kimb", tau[a], dP[..., b], optimize=True)
            )
    return P, dP, d2P


def sigma_from_geometry(H_of, P, dP, d2P, kpts, omegas, components, nocc):
    """sigma from the projector C-trace (arXiv:2409.16358 Eq. (6) family),
    with the SAME assembly as the Stage 26 oracle (C_sym, (bc)-sym, f_nm delta,
    -2pi/Nk Im)."""
    nk = len(kpts)
    energy, _ = np.linalg.eigh(H_of(kpts))[:2]
    C = np.zeros((nk, NW, NW, 3, 3, 3), dtype=np.complex128)
    for m in range(NW):
        for n in range(NW):
            Pn = P[..., n]
            for a in range(3):
                for c in range(3):
                    Z1 = d2P[..., n, a, c]
                    Z2 = np.einsum("kij,kjm->kim", dP[..., m, a], dP[..., n, c], optimize=True)
                    W = Z1 + Z2
                    for b in range(3):
                        Y = dP[..., m, b]
                        C[:, m, n, a, b, c] = np.einsum(
                            "kij,kjl,kli->k", Pn, Y, W, optimize=True
                        )
    C_sym = 0.5 * (C + np.swapaxes(C, 4, 5))
    occupation = np.zeros(NW, dtype=np.float64)
    occupation[:nocc] = 1.0
    sigma = np.zeros((len(components), len(omegas)), dtype=np.float64)
    for ci, (a1, b1, c1) in enumerate(components):
        a, b, c = a1 - 1, b1 - 1, c1 - 1
        for wi, omega in enumerate(omegas):
            for n in range(NW):
                for m in range(NW):
                    if m <= n:
                        continue
                    f_nm = occupation[n] - occupation[m]
                    eps_mn = energy[:, m] - energy[:, n]
                    delta = ETA_SP / np.pi / ((omega - eps_mn) ** 2 + ETA_SP ** 2)
                    sigma[ci, wi] += float(np.sum(
                        f_nm * delta * (-2.0 * np.pi / nk)
                        * np.imag(C_sym[:, n, m, a, b, c])
                    ))
    return sigma


def raw_from_band_matrices(energy, vmn, wmnab, nocc, eta_band, omegas, components):
    """Exact mirror of NLO `_compute_block` (identical to oracle model_raw tail)."""
    nk = vmn.shape[0]
    omn = energy[:, :, None] - energy[:, None, :]
    factor = -1j * omn / (omn * omn + eta_band * eta_band)
    diag = np.arange(NW)
    factor[:, diag, diag] = 0
    rmn = vmn * factor[..., None]
    eta2 = eta_band * eta_band
    rmndk = np.zeros((nk, NW, NW, 3, 3), dtype=np.complex128)
    for m in range(NW):
        for n in range(NW):
            if m == n:
                continue
            dv_ab = np.zeros((nk, 3, 3), dtype=np.complex128)
            int_ab = np.zeros((nk, 3, 3), dtype=np.complex128)
            for p in range(NW):
                if p != m and p != n:
                    int_ab += (
                        vmn[:, n, p, :, None] * vmn[:, p, m, None, :]
                        * (omn[:, p, m] / (omn[:, p, m] ** 2 + eta2))[:, None, None]
                        - vmn[:, n, p, :, None] * vmn[:, p, m, None, :]
                        * (omn[:, n, p] / (omn[:, n, p] ** 2 + eta2))[:, None, None]
                    )
            for a in range(3):
                for b in range(3):
                    dv_ab[:, a, b] = (
                        vmn[:, n, m, a] * (vmn[:, n, n, b] - vmn[:, m, m, b])
                        + vmn[:, n, m, b] * (vmn[:, n, n, a] - vmn[:, m, m, a])
                    )
            gap = omn[:, n, m]
            den = gap * gap + eta2
            rmndk[:, n, m] = (1j * gap / den)[:, None, None] * (
                dv_ab * gap[:, None, None] / den[:, None, None]
                - wmnab[:, n, m]
                + int_ab
            )
    occupation = np.zeros(NW, dtype=np.float64)
    occupation[:nocc] = 1.0
    occdiff = occupation[None, :] - occupation[:, None]
    kweight = 1.0 / nk
    raw = np.zeros((len(components), len(omegas)), dtype=np.float64)
    for ci, (a1, b1, c1) in enumerate(components):
        a, b, c = a1 - 1, b1 - 1, c1 - 1
        iabc = rmn[..., b] * np.swapaxes(rmndk[..., c, a], 1, 2)
        iacb = rmn[..., c] * np.swapaxes(rmndk[..., b, a], 1, 2)
        weights = occdiff * (iabc + iacb)
        for wi, omega in enumerate(omegas):
            delta = ETA_SP / np.pi / ((omn - omega) ** 2 + ETA_SP ** 2)
            raw[ci, wi] = float(np.real(1j * kweight * np.sum(weights * delta)))
    return raw


def raw_source_gauge(kpts, nocc, eta_band, omegas, components):
    """R_src: raw from SOURCE-gauge analytic H_src and its derivatives."""
    nk = len(kpts)
    energy = np.zeros((nk, NW), dtype=np.float64)
    vmn = np.zeros((nk, NW, NW, 3), dtype=np.complex128)
    wmnab = np.zeros((nk, NW, NW, 3, 3), dtype=np.complex128)
    dH, d2H = d1_d2_src_at(kpts)
    for i in range(nk):
        H = h_src_at(kpts)[i]
        ev, wf = np.linalg.eigh(H)
        energy[i] = ev
        vmn[i] = np.einsum("in,ijd,jm->nmd", wf.conj(), dH[i], wf, optimize=True)
        t = np.einsum("in,ijab,jm->nmab", wf.conj(), d2H[i], wf, optimize=True)
        wmnab[i] = t
    return raw_from_band_matrices(energy, vmn, wmnab, nocc, eta_band, omegas, components)


def dev_stats(ref, val):
    dev = np.abs(ref - val)
    strong = np.abs(ref) > STRONG_RAW
    return {
        "max_abs_dev_all": float(dev.max()),
        "max_abs_dev_strong": float(dev[strong].max()) if strong.any() else None,
        "n_strong_points": int(strong.sum()),
    }


def mean_ratio(ref, val):
    """mean val/ref over strong points (reports the measured constant)."""
    sel = np.abs(ref) > STRONG_RAW
    if not sel.any():
        return None
    return float(np.mean((val / ref)[sel]))


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else _HERE.parent
    root = Path(root).resolve()
    # optional fixture override for third-party reproduction:
    #   python3 stage27_gauge_covariance_audit.py <output-root> <fixture-dir>
    if len(sys.argv) > 2:
        fixture_dir = Path(sys.argv[2]).resolve()
    else:
        fixture_dir = MODEL_DIR
    centers = read_centers(fixture_dir / "wannier90-wcenter.dat")
    blocks = read_hr(fixture_dir / "wannier90_hr.dat")

    kpts5 = grid_kpoints((5, 5, 1))

    # ---- G0: analytic H_cell == HR kernel H(k) ----
    Hk_hr = fourier_for_sigma(blocks, None, kpts5)
    g0 = float(np.max(np.abs(Hk_hr - h_cell_at(kpts5))))

    # ---- G0b: tau from U == center differences from the package ----
    tau = tau_from_U(kpts5)
    delta_c = centers[1] - centers[0]  # (1.0, 0.0, 0.0) A
    tau_diag = np.array([tau[a][1, 1] for a in range(3)])
    offdiag = np.array(
        [np.abs(tau[a] - np.diag(np.diag(tau[a]))).max() for a in range(3)]
    )
    g0b = {
        "tau_from_U_diag": [float(np.real(tau[a][1, 1])) for a in range(3)],
        "tau_from_U_imag_diag": [float(np.imag(tau[a][1, 1])) for a in range(3)],
        "center_diff_A": [float(x) for x in delta_c],
        "offdiag_abs_max": float(offdiag.max()),
        "max_dev_diag": float(np.max(np.abs(tau_diag - delta_c))),
    }

    # ---- G1: D_a H_cell == U^dagger (d_a H_src) U ----
    U0 = u_at(kpts5)
    Ud0 = U0.conj().transpose(0, 2, 1)
    dH_src, _ = d1_d2_src_at(kpts5)
    rot_dH = np.einsum("kij,kjmd,kmn->kind", Ud0, dH_src, U0, optimize=True)
    g1 = 0.0
    for a in range(3):
        sa = reduced_shift(H_D1 * np.eye(3)[a])
        dHc = (h_cell_at(kpts5 + sa) - h_cell_at(kpts5 - sa)) / (2 * H_D1)
        Hc0 = h_cell_at(kpts5)
        D_Hc = dHc + 1j * (np.einsum("kij,jm->kim", Hc0, tau[a], optimize=True)
                           - np.einsum("ij,kjm->kim", tau[a], Hc0, optimize=True))
        g1 = max(g1, float(np.max(np.abs(D_Hc - rot_dH[..., a]))))

    # ---- G2: D_a P_cell == U^dagger (d_a P_src) U ----
    def P_src_at(ks):
        return projectors_from_H(h_src_at, ks)

    P_src0 = P_src_at(kpts5)
    rot_dP_all = np.zeros((len(kpts5), NW, NW, NW, 3), dtype=np.complex128)
    for a in range(3):
        sa = reduced_shift(H_D1 * np.eye(3)[a])
        dPs = (P_src_at(kpts5 + sa) - P_src_at(kpts5 - sa)) / (2 * H_D1)
        rot_dP_all[..., a] = np.einsum("kij,kjmb,kom->kiob", Ud0, dPs, U0, optimize=True)
    g2 = 0.0
    P_cell0 = projectors_from_H(h_cell_at, kpts5)
    for a in range(3):
        sa = reduced_shift(H_D1 * np.eye(3)[a])
        dPc = (projectors_from_H(h_cell_at, kpts5 + sa)
               - projectors_from_H(h_cell_at, kpts5 - sa)) / (2 * H_D1)
        D_Pc = dPc + 1j * (np.einsum("kijb,jm->kimb", P_cell0, tau[a], optimize=True)
                           - np.einsum("ij,kjmb->kimb", tau[a], P_cell0, optimize=True))
        g2 = max(g2, float(np.max(np.abs(D_Pc - rot_dP_all[..., a]))))

    # ---- G3/G4: sigma geometries (identical assembly, unified conventions) ----
    H_src_of = lambda ks: h_src_at(ks)
    H_cell_of = lambda ks: h_cell_at(ks)

    P_s, dP_s, d2P_s = plain_geometry(H_src_of, kpts5)
    sig_src = sigma_from_geometry(H_src_of, P_s, dP_s, d2P_s, kpts5, OMEGAS, COMPONENTS, 1)

    P_c, dP_c, d2P_c = plain_geometry(H_cell_of, kpts5)
    sig_cell_plain = sigma_from_geometry(H_cell_of, P_c, dP_c, d2P_c, kpts5, OMEGAS, COMPONENTS, 1)

    tau = tau_from_U(kpts5)
    P_cov, dP_cov, d2P_cov = covariant_geometry(H_cell_of, kpts5, tau)
    sig_cell_cov = sigma_from_geometry(H_cell_of, P_cov, dP_cov, d2P_cov, kpts5, OMEGAS, COMPONENTS, 1)

    # ---- raws (same omegas/components/eta_sp/eta_band/k mesh) ----
    R_src = raw_source_gauge(kpts5, 1, ETA_BAND_PRIMARY, OMEGAS, COMPONENTS)
    R_cell_plain = model_raw(blocks, kpts5, 1, ETA_BAND_PRIMARY, OMEGAS, COMPONENTS,
                             centers=None, legacy=False)
    R_dressed = model_raw(blocks, kpts5, 1, ETA_BAND_PRIMARY, OMEGAS, COMPONENTS,
                          centers=centers, legacy=False)

    checks = {
        "G0_Hcell_analytic_vs_HR_kernel": g0,
        "G0b_tau_vs_wcenter": g0b,
        "G1_DHcell_eq_Udag_dHsrc_U": g1,
        "G2_DPcell_eq_Udag_dPsrc_U": g2,
        "G3_sigma_src_vs_cell_cov": dev_stats(sig_src, sig_cell_cov),
        "G4_sigma_cell_plain_vs_cell_cov": dev_stats(sig_cell_cov, sig_cell_plain),
        "G5_sigma_cell_cov_vs_minus_pi_Rdressed": dev_stats(-np.pi * R_dressed, sig_cell_cov),
        "G6_sigma_cell_plain_vs_minus_pi_Rplain": dev_stats(-np.pi * R_cell_plain, sig_cell_plain),
        "G7_Rdressed_vs_Rsrc": float(np.max(np.abs(R_dressed - R_src))),
        "G8_sigma_src_vs_minus_pi_Rsrc": dev_stats(-np.pi * R_src, sig_src),
    }

    # ---- pairing table (honest documentation of ALL pairings) ----
    sigmas = {"src": sig_src, "cell_plain": sig_cell_plain, "cell_cov": sig_cell_cov}
    raws = {"Rsrc": R_src, "Rplain": R_cell_plain, "Rdressed": R_dressed}
    pairing = {}
    for sname, s in sigmas.items():
        for rname, r in raws.items():
            pairing[f"{sname}_over_pi_{rname}"] = {
                "max_abs_dev_strong": dev_stats(-np.pi * r, s)["max_abs_dev_strong"],
                "mean_ratio_sigma_over_minus_pi_raw": mean_ratio(-np.pi * r, s),
            }

    # ---- G5b: eta_band stability and k-grid convergence ----
    eta_series = {}
    for eta in ETA_BAND_GRID:
        R = model_raw(blocks, kpts5, 1, eta, OMEGAS, COMPONENTS, centers=centers, legacy=False)
        eta_series[str(eta)] = dev_stats(-np.pi * R, sig_cell_cov)
    kgrid_series = {}
    for n in KGRID:
        kp = grid_kpoints((n, n, 1))
        R = model_raw(blocks, kp, 1, ETA_BAND_PRIMARY, OMEGAS, COMPONENTS, centers=centers,
                      legacy=False)
        Pn, dPn, d2Pn = covariant_geometry(H_cell_of, kp, tau_from_U(kp))
        sg = sigma_from_geometry(H_cell_of, Pn, dPn, d2Pn, kp, OMEGAS, COMPONENTS, 1)
        kgrid_series[str(n)] = dev_stats(-np.pi * R, sg)

    result = {
        "schema": "stage27-gauge-covariance-audit/v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "conventions": {
            "H_src": "analytic 2-band Haldane M=0.4 phi=pi/2 t1=1.0 t2=0.2 NN=(3,2) NNN=(3,2)",
            "U": "diag(1, e^{-ik.d1}), d1=(1,0) A (Cartesian-k convention)",
            "tau": "i(d_a U)U^dagger, verified against wannier90-wcenter.dat",
            "covariant_derivative": "D_a X = d_a X + i[X, tau_a]",
            "k_grid": "5x5x1 (primary), 9x9x1/13x13x1 (convergence)",
            "eta_sp": ETA_SP,
            "eta_band_primary": ETA_BAND_PRIMARY,
            "eta_band_scan": list(ETA_BAND_GRID),
            "components": list(COMPONENTS),
            "broadening_delta": "eta_sp/pi/((w-eps)^2+eta_sp^2)",
            "bc_sym": "C_sym = (C + swapaxes(b,c))/2, sigma from Im C_sym",
            "projector_fd_steps": {"h_d1": H_D1, "h_d2": H_D2},
            "strong_point_mask": f"|raw| > {STRONG_RAW}",
        },
        "checks": checks,
        "pairing_table": pairing,
        "eta_band_stability_G5": eta_series,
        "kgrid_convergence_G5": kgrid_series,
        "measured_constants": {
            "sigma_cell_cov_over_minus_Rdressed": mean_ratio(-R_dressed, sig_cell_cov),
            "sigma_src_over_minus_Rsrc": mean_ratio(-R_src, sig_src),
        },
        "negative_findings_kept": [
            "G4: plain cell-gauge sigma differs from the covariant/source sigma "
            "(ordinary dP_cell is NOT the physically matched object); the plain "
            "pair is gauge-dependent.",
            "Stage 25's -(pi/2) relation is NOT assumed; the measured constant is "
            "~pi for the matched (covariant/dressed, source) pairings, see pairing table.",
            "Pairings that mix conventions (e.g. cell_plain sigma with dressed raw) "
            "are reported in the pairing table, not hidden.",
        ],
        "unresolved_items": [
            "Absolute SI prefactor and 2D sheet normalization (Lz*0.1) remain out of "
            "scope of this geometric audit.",
            "FDR SM Eq. (75) original text not available on this machine "
            "(SKIPPED_NO_FDR_ARCHIVE); the projector form is an independent "
            "reconstruction of arXiv:2409.16358 Eqs. (3)+(6), not a verbatim SM copy.",
        ],
    }
    out = root / "STAGE27_GAUGE_COVARIANCE_AUDIT.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in ("pairing_table",)},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
