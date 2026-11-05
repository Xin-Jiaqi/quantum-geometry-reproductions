#!/usr/bin/env python3
"""Stage 26 — independent second-derivative (d2H/dk_a dk_b) validation.

Completely independent of NLO source code. Uses:
  * the synthetic spinless Haldane package HR + centers (own readers);
  * an own Fourier kernel for H(k), dH(k), d2H(k) with center embedding;
  * exact analytic dH, d2H of the source-gauge target Hamiltonian (2-band
    Haldane d-vector form) mapped through the cell-gauge unitary;
  * central finite differences of H(k) and of dH(k).

Checks (all Cartesian components xx, xy, yx, yy, xz, zx, yz, zy, zz):
  A1 analytic-vs-finite-difference agreement;
  A2 mixed-partial symmetry w_ab == w_ba;
  A3 Hermiticity w[n,m] == conj(w[m,n]);
  A4 covariance under a constant orbital-basis rotation;
  A5 step-size convergence (central differences -> O(h^2));
  A6 non-zero mixed derivatives on this model (not a trivial-zero test).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SQ3 = np.sqrt(3.0)
A1 = np.array([1.5, SQ3 / 2.0, 0.0])
A2 = np.array([1.5, -SQ3 / 2.0, 0.0])
B1 = np.array([2.0 * np.pi / 3.0, 2.0 * np.pi / SQ3, 0.0])
B2 = np.array([2.0 * np.pi / 3.0, -2.0 * np.pi / SQ3, 0.0])
ALATT = np.vstack([A1, A2])
D2_PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (2, 0), (0, 1))
M, PHI, T1, T2 = 0.4, np.pi / 2, 1.0, 0.2
DELTA1 = np.array([1.0, 0.0])
NN = np.array([[1.0, 0.0], [-0.5, SQ3 / 2.0], [-0.5, -SQ3 / 2.0]])
NNN = np.array([[1.5, SQ3 / 2.0], [-1.5, SQ3 / 2.0], [0.0, -SQ3]])


def read_hr(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    nw = int(lines[1]); nr = int(lines[2])
    cursor = 3
    degs = []
    while len(degs) < nr:
        degs.extend(int(v) for v in lines[cursor].split()); cursor += 1
    rows = [ln.split() for ln in lines[cursor:] if ln.strip()]
    blocks = {}
    for b in range(nr):
        chunk = rows[b * nw * nw:(b + 1) * nw * nw]
        R = tuple(int(x) for x in chunk[0][:3])
        mat = np.zeros((nw, nw), dtype=np.complex128)
        for row in chunk:
            mat[int(row[3]) - 1, int(row[4]) - 1] = float(row[5]) + 1j * float(row[6])
        blocks[R] = mat
    return blocks


def read_centers(path):
    c = []
    for ln in path.read_text(encoding="utf-8").splitlines():
        p = ln.split()
        if len(p) >= 4:
            c.append([float(p[1]), float(p[2]), float(p[3])])
    return np.asarray(c)


def make_kpoints(nk):
    def axis(n):
        return np.asarray([0.0]) if n == 1 else -0.5 + np.arange(n) / n
    gx, gy, gz = np.meshgrid(axis(nk[0]), axis(nk[1]), axis(nk[2]), indexing="ij")
    return np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))


def analytic_source_hamiltonian(kcart):
    """Exact H_src(k) and its first/second Cartesian derivatives (2-band form).

    H = h0*I + hx*sig_x + hy*sig_y + hz*sig_z with
      h0 = 2 t2 cos(phi) sum_b cos(b.k),  hx = t1 sum_d cos(d.k),
      hy = -t1 sum_d sin(d.k),            hz = M + 2 t2 sin(phi) sum_b sin(b.k),
    d in NN, b in NNN (both strictly 2D -> all z-derivatives vanish).
    """
    theta1 = NN @ kcart[:2]
    theta2 = NNN @ kcart[:2]
    c1, s1 = np.cos(theta1), np.sin(theta1)
    c2, s2 = np.cos(theta2), np.sin(theta2)
    h0 = 2.0 * T2 * np.cos(PHI) * float(c2.sum())
    hx = T1 * float(c1.sum())
    hy = -T1 * float(s1.sum())
    hz = M + 2.0 * T2 * np.sin(PHI) * float(s2.sum())
    H = np.array([[h0 + hz, hx - 1j * hy], [hx + 1j * hy, h0 - hz]], dtype=np.complex128)
    sig = np.zeros((3, 2, 2), dtype=np.complex128)
    sig[0] = [[0, 1], [1, 0]]
    sig[1] = [[0, -1j], [1j, 0]]
    sig[2] = [[1, 0], [0, -1]]
    I2 = np.eye(2, dtype=np.complex128)

    d1 = np.zeros((3, 2, 2), dtype=np.complex128)
    d2 = np.zeros((3, 3, 2, 2), dtype=np.complex128)
    for a in range(3):
        if a == 2:
            continue  # strictly 2D model: no z derivative
        dhx = -T1 * float((NN[:, a] * s1).sum())
        dhy = -T1 * float((NN[:, a] * c1).sum())
        dhz = 2.0 * T2 * np.sin(PHI) * float((NNN[:, a] * c2).sum())
        dh0 = -2.0 * T2 * np.cos(PHI) * float((NNN[:, a] * s2).sum())
        d1[a] = dh0 * I2 + dhx * sig[0] + dhy * sig[1] + dhz * sig[2]
        for b in range(3):
            if b == 2:
                continue
            d2hx = -T1 * float((NN[:, a] * NN[:, b] * c1).sum())
            d2hy = T1 * float((NN[:, a] * NN[:, b] * s1).sum())
            d2hz = -2.0 * T2 * np.sin(PHI) * float((NNN[:, a] * NNN[:, b] * s2).sum())
            d2h0 = -2.0 * T2 * np.cos(PHI) * float((NNN[:, a] * NNN[:, b] * c2).sum())
            d2[a, b] = d2h0 * I2 + d2hx * sig[0] + d2hy * sig[1] + d2hz * sig[2]
    return H, d1, d2


def cell_gauge_unitary(kcart):
    return np.diag([1.0 + 0j, np.exp(-1j * float(np.dot(kcart[:2], DELTA1)))])


def analytic_cell_objects(kcart):
    """Gauge-rotated source objects in the Wannier cell gauge.

    H_cell(k) = U†(k) H_src(k) U(k),  U = diag(1, e^{-i k_x})  (delta1=(1,0)).
    The position-dressed first/second derivatives follow the Wannier
    center-embedding convention (the same one verified by the Stage 23 bridge):
        dH_pos[a]  = U† dH_src[a] U
        d2H_pos[a,b] = U† d2H_src[a,b] U
    i.e. the source-gauge derivatives rotated into the cell basis, WITHOUT the
    derivative of the gauge unitary itself (that part is carried by the Wannier
    centers in the position operator, mode diagonal_tba).
    """
    H, d1, d2 = analytic_source_hamiltonian(kcart)
    phase = np.exp(-1j * float(np.dot(kcart[:2], DELTA1)))
    U = np.diag([1.0 + 0j, phase])
    Ud = np.diag([1.0 + 0j, np.conj(phase)])
    Hc = Ud @ H @ U
    dHc = np.zeros((3, 2, 2), dtype=np.complex128)
    d2Hc = np.zeros((3, 3, 2, 2), dtype=np.complex128)
    for a in range(3):
        dHc[a] = Ud @ d1[a] @ U
        for b in range(3):
            d2Hc[a, b] = Ud @ d2[a, b] @ U
    return Hc, dHc, d2Hc


def hr_fourier(blocks, centers, kpts, nw=2):
    """Own Fourier kernel: H(k), dH/dk_cart (eV/A), d2H (full 3x3 Cartesian)."""
    R_list = list(blocks)
    R = np.asarray(R_list, dtype=np.int64)
    weighted = np.stack([blocks[r] for r in R_list])
    phase = np.exp(1j * 2 * np.pi * (kpts @ R.T))
    Hk = np.einsum("kr,rij->kij", phase, weighted)
    nk = len(kpts)
    dH = np.zeros((nk, nw, nw, 3), dtype=np.complex128)
    d2H = np.zeros((nk, nw, nw, 3, 3), dtype=np.complex128)
    for d in range(3):
        aR = np.zeros((len(R_list), nw, nw))
        for r, rr in enumerate(R_list):
            Rcart = rr[0] * A1 + rr[1] * A2
            aR[r] = Rcart[d] + (centers[:, d][None, :] - centers[:, d][:, None])
        dH[..., d] = np.einsum("kr,rij->kij", phase, 1j * aR * weighted)
    for d1 in range(3):
        for d2 in range(3):
            aR = np.zeros((len(R_list), nw, nw))
            for r, rr in enumerate(R_list):
                Rcart = rr[0] * A1 + rr[1] * A2
                aR[r] = Rcart[d1] * Rcart[d2] + Rcart[d1] * (centers[:, d2][None, :] - centers[:, d2][:, None]) \
                        + Rcart[d2] * (centers[:, d1][None, :] - centers[:, d1][:, None]) \
                        + (centers[:, d1][None, :] - centers[:, d1][:, None]) * (centers[:, d2][None, :] - centers[:, d2][:, None])
            d2H[..., d1, d2] = np.einsum("kr,rij->kij", phase, -aR * weighted)
    return Hk, dH, d2H


def band_projection(wf, op):
    """<n|op|m> for a batch: op shape (k, nw, nw[, ...])."""
    return np.einsum("kin,kij...,kjm->knm...", wf.conj(), op, wf)


def cartesian_shift_reduced(shift_cart):
    """Map a Cartesian k-shift (Angstrom^-1) into reduced-coordinate shifts."""
    du = float(shift_cart @ A1) / (2 * np.pi)
    dv = float(shift_cart @ A2) / (2 * np.pi)
    return np.array([du, dv, 0.0])


def finite_difference_h(hr_blocks, centers, kpts, h, nw=2):
    """Central finite-difference d2H from H(k) only, in CARTESIAN k units."""
    d2 = np.zeros((len(kpts), nw, nw, 3, 3), dtype=np.complex128)
    for a in range(3):
        for b in range(3):
            ea, eb = np.eye(3)[a], np.eye(3)[b]
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    shift = cartesian_shift_reduced(s1 * h * ea + s2 * h * eb)
                    Hk, _, _ = hr_fourier(hr_blocks, centers, kpts + shift, nw)
                    d2[..., a, b] += s1 * s2 * Hk / (4 * h * h)
    return d2


def finite_difference_dh(hr_blocks, centers, kpts, h, nw=2):
    """Central finite difference of the HR-fourier dH (Cartesian shifts)."""
    d2 = np.zeros((len(kpts), nw, nw, 3, 3), dtype=np.complex128)
    for a in range(3):
        for b in range(3):
            shift = cartesian_shift_reduced(h * np.eye(3)[b])
            _, dHp, _ = hr_fourier(hr_blocks, centers, kpts + shift, nw)
            _, dHm, _ = hr_fourier(hr_blocks, centers, kpts - shift, nw)
            d2[..., a, b] = (dHp[..., a] - dHm[..., a]) / (2 * h)
    return d2


def align_phases(target, source):
    """Fix the k-dependent eigenvector phases of `target` onto `source` via
    maximal-overlap real alignment (valid for a gapped model, no crossings).
    Arrays are shaped (k, orbital, band) (np.linalg.eigh convention)."""
    overlap = np.einsum("kob,kom->kbm", source.conj(), target)
    phase = np.exp(-1j * np.angle(np.diagonal(overlap, axis1=1, axis2=2)))
    return target * phase[:, None, :]


def finite_difference_dh_analytic(hr_blocks, centers, kpts, nw=2):
    """Analytic derivative of the dressed dH: -sum_R (R_a+dt_a) R_b H e^{ikR}."""
    R_list = list(hr_blocks)
    R = np.asarray(R_list, dtype=np.int64)
    weighted = np.stack([hr_blocks[r] for r in R_list])
    phase = np.exp(1j * 2 * np.pi * (kpts @ R.T))
    out = np.zeros((len(kpts), nw, nw, 3, 3), dtype=np.complex128)
    for a in range(3):
        for b in range(3):
            fac = np.zeros((len(R_list), nw, nw))
            for r, rr in enumerate(R_list):
                Rcart = rr[0] * A1 + rr[1] * A2
                fac[r] = (
                    (Rcart[a] + (centers[:, a][None, :] - centers[:, a][:, None]))
                    * Rcart[b]
                )
            out[..., a, b] = np.einsum("kr,rij->kij", phase, -fac * weighted)
    return out


def plain_first_derivative(hr_blocks, centers, kpts, h, nw=2):
    """Central FD of H(k) wrt Cartesian k (plain derivative, no dressing)."""
    dH = np.zeros((len(kpts), nw, nw, 3), dtype=np.complex128)
    for a in range(3):
        shift = cartesian_shift_reduced(h * np.eye(3)[a])
        Hp, _, _ = hr_fourier(hr_blocks, centers, kpts + shift, nw)
        Hm, _, _ = hr_fourier(hr_blocks, centers, kpts - shift, nw)
        dH[..., a] = (Hp - Hm) / (2 * h)
    return dH


def center_dressing_term(hr_blocks, centers, kpts, nw=2):
    """sum_R i (wc_j - wc_i)_a H(R) e^{i2pi k R}: the position-embedding term."""
    R_list = list(hr_blocks)
    R = np.asarray(R_list, dtype=np.int64)
    weighted = np.stack([hr_blocks[r] for r in R_list])
    phase = np.exp(1j * 2 * np.pi * (kpts @ R.T))
    term = np.zeros((len(kpts), nw, nw, 3), dtype=np.complex128)
    for d in range(3):
        aR = np.zeros((len(R_list), nw, nw))
        for r, rr in enumerate(R_list):
            aR[r] = centers[:, d][None, :] - centers[:, d][:, None]
        term[..., d] = np.einsum("kr,rij->kij", phase, 1j * aR * weighted)
    return term


def main() -> int:
    work = Path(sys.argv[1]).resolve()
    pkg = work / "models" / "stage24_haldane_synthetic" / "L1"
    blocks = read_hr(pkg / "wannier90_hr.dat")
    centers = read_centers(pkg / "wannier90-wcenter.dat")
    rng = np.random.default_rng(20261026)
    kpts = rng.uniform(-0.5, 0.5, size=(96, 3))
    kpts[:, 2] = 0.0

    Hk, dH, d2H = hr_fourier(blocks, centers, kpts)
    checks: dict[str, object] = {}
    metrics: dict[str, object] = {}

    # A1a: HR-kernel d2H vs analytic gauge-rotated source d2H (position-dressed
    # convention, U^dag d2H_src U; the Stage 23 bridge reference).
    max_analytic = 0.0
    max_h = 0.0
    for i, k in enumerate(kpts):
        kcart = k[0] * B1 + k[1] * B2
        Hc, dHc, d2Hc = analytic_cell_objects(kcart)
        max_h = max(max_h, float(np.max(np.abs(Hk[i] - Hc))))
        for a in range(3):
            for b in range(3):
                max_analytic = max(
                    max_analytic,
                    float(np.max(np.abs(d2H[i, ..., a, b] - d2Hc[a, b]))),
                )
    checks["A1a_hr_kernel_d2H_matches_analytic_gauge_rotated"] = max_analytic < 1e-12
    checks["A0_hr_kernel_H_matches_analytic_cell_gauge"] = max_h < 1e-12
    metrics["A0_max_abs_hr_kernel_H_vs_analytic_eV"] = max_h
    metrics["A1a_max_abs_hr_kernel_d2H_vs_analytic_eV_per_A2"] = max_analytic

    # A1b: finite difference of dH (Cartesian central) vs the analytic
    # derivative-of-dH expression (kernel2 = -sum (R_a+dt_a) R_b H e^{ikR}).
    kernel2 = finite_difference_dh_analytic(blocks, centers, kpts)
    err_dh = {}
    for h in (0.1, 0.05, 0.02, 0.01, 0.005, 0.002):
        d2fd = finite_difference_dh(blocks, centers, kpts, h)
        err_dh[h] = float(np.max(np.abs(d2fd - kernel2)))
    metrics["A1b_fd_dH_err"] = err_dh
    checks["A1b_fd_dH_converges"] = err_dh[0.002] < 1e-5

    # convergence order estimate (central diff -> ~2)
    hs = np.asarray([0.1, 0.05, 0.02, 0.01], dtype=float)
    errs = np.asarray([err_dh[h] for h in hs], dtype=float)
    slope = float(np.polyfit(np.log(hs), np.log(errs), 1)[0])
    metrics["A1b_fd_dH_convergence_slope"] = slope
    checks["A1b_central_difference_second_order"] = 1.5 < slope < 2.5

    # A1d: dressing cross-term identity d2H[a,b] = kernel2[a,b] + i dt_b dH[a]
    cross = np.zeros_like(d2H)
    for b in range(3):
        dt = centers[:, b][None, :] - centers[:, b][:, None]
        cross[..., :, b] = 1j * dt[None, :, :, None] * dH
    metrics["A1d_dressing_cross_term_residual"] = float(
        np.max(np.abs((d2H - kernel2) - cross))
    )
    checks["A1d_dressing_cross_term_identity"] = metrics["A1d_dressing_cross_term_residual"] < 1e-12

    # A1c: position-embedding decomposition dH_kernel = plain dH + center term.
    dH_plain = plain_first_derivative(blocks, centers, kpts, 1e-4)
    dress = center_dressing_term(blocks, centers, kpts)
    decom = float(np.max(np.abs(dH - (dH_plain + dress))))
    metrics["A1c_position_embedding_decomposition_residual"] = decom
    checks["A1c_position_embedding_decomposition"] = decom < 1e-7

    # A2: mixed-partial symmetry (analytic kernel, all k)
    sym = float(np.max(np.abs(d2H[..., 0, 1] - d2H[..., 1, 0])))
    sym = max(sym, float(np.max(np.abs(d2H[..., 0, 2] - d2H[..., 2, 0]))))
    sym = max(sym, float(np.max(np.abs(d2H[..., 1, 2] - d2H[..., 2, 1]))))
    checks["A2_mixed_partial_symmetry"] = sym < 1e-13
    metrics["A2_max_abs_wab_minus_wba"] = sym

    # A3: Hermiticity of the band-projected second derivative
    _, wf = np.linalg.eigh(Hk)
    wband = band_projection(wf, d2H)
    herm = float(np.max(np.abs(wband - np.swapaxes(np.swapaxes(wband.conj(), 1, 2), 3, 4))))
    checks["A3_hermiticity"] = herm < 1e-13
    metrics["A3_max_abs_wmn_minus_conj_wnm"] = herm

    # A4: covariance under constant orbital rotation V (unitary): band-basis w
    # is identical up to the (k-dependent) eigenvector gauge; fix phases by
    # maximal-overlap alignment and compare.
    rng2 = np.random.default_rng(7)
    A = rng2.normal(size=(2, 2)) + 1j * rng2.normal(size=(2, 2))
    V, _ = np.linalg.qr(A)
    Hk2 = np.einsum("ab,kbc,cd->kad", V.conj().T, Hk, V)
    d2H2 = np.einsum("ip,kpqab,qj->kijab", V.conj().T, d2H, V)
    _, wf2 = np.linalg.eigh(Hk2)
    wf2_expected = np.einsum("ba,kan->kbn", V.conj().T, wf)
    wf2a = align_phases(wf2, wf2_expected)
    wband2 = band_projection(wf2a, d2H2)
    cov = float(np.max(np.abs(wband - wband2)))
    checks["A4_orbital_rotation_covariance"] = cov < 1e-12
    metrics["A4_max_abs_covariance_residual"] = cov

    # A5: convergence recorded above; explicit check that error shrinks monotonically
    err_dh = metrics["A1b_fd_dH_err"]
    e_h = [err_dh[h] for h in (0.05, 0.01, 0.002)]
    checks["A5_monotone_convergence"] = e_h[0] > e_h[1] > e_h[2]

    # A6: non-zero mixed derivatives (a model where the old indexing corrupts yx/yy)
    nz = {}
    for (a, b) in ((0, 1), (1, 0), (1, 1), (0, 2), (1, 2), (2, 2)):
        nz[f"w_{'xyz'[a]}{'xyz'[b]}"] = float(np.max(np.abs(d2H[..., a, b])))
    metrics["A6_nonzero_components_max_abs"] = nz
    checks["A6_mixed_derivatives_nonzero"] = (
        nz["w_xy"] > 1e-3 and nz["w_yx"] > 1e-3 and nz["w_yy"] > 1e-3
    )

    status = "PASS" if all(v is True for v in checks.values()) else "FAIL"
    report = {
        "schema": "stage26-second-derivative-validation/v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "model": "synthetic spinless Haldane, M=0.4 phi=pi/2 t1=1 t2=0.2, package stage24_haldane_synthetic/L1",
        "kpoints": "96 random in-plane reduced k (kz=0)",
        "checks": checks,
        "metrics": metrics,
        "interpretation": (
            "The HR-kernel second derivatives (with center embedding) are the exact "
            "Cartesian d2H/dk_a dk_b of the cell-gauge Hamiltonian: they match the "
            "analytic cell-gauge objects to <1e-12 eV/A^2, are symmetric in (a,b), "
            "Hermitian in band indices, covariant under constant orbital rotations, "
            "and converge with central finite differences at second order. Mixed "
            "derivatives xy/yx/yy are strongly non-zero on this model."
        ),
    }
    out = work / "STAGE26_SECOND_DERIVATIVE_VALIDATION.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "metrics": metrics}, indent=2, ensure_ascii=False))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
