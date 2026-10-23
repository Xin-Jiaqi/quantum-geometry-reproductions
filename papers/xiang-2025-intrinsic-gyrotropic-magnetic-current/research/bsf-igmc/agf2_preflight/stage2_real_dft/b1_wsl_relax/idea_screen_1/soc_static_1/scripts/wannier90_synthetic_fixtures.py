#!/usr/bin/env python3
"""Deterministic SYNTHETIC authoring of genuine Wannier90 file syntax.

NO Wannier90 run. Labelled synthetic-only fixtures that validate the strict
parser (import_hsa_operator_bundle.py) against known operators. Layouts
follow the official Wannier90 documentation (master, fetched 2026-10-23)
and the testw90_rmn benchmark r.dat:

  hr.dat : comment; num_wann; nrpts; degeneracy block (15/line); per
           (R, m, n): R1 R2 R3 m n Re Im   (m fastest, 1-based)
  r.dat  : comment; num_wann; nrpts; per (R, m, n):
           R1 R2 R3 m n xRe xIm yRe yIm zRe zIm (no degeneracy block)
  spn    : comment; num_bands num_kpts; per k: upper triangle
           <psi_n|sigma_i|psi_m>, m=1..nb, n=1..m, i=x,y,z (i fastest)
  u.mat  : comment; num_kpts num_wann num_wann; blank; per k: kpoint line +
           num_wann^2 pairs, column-major, blank between
  u_dis.mat / v.mat : comment; num_kpts num_wann num_bands; same layout
  eig    : 'band k_ordinal eigenvalue' (1-based)
  win    : explicit kpoints block, saxis, spinors, unit_cell_cart
           (orthorhombic 3x3), write_ndegen_applied = true

Gauge algebra (official W90, contract doc §2):
  V(k) = U_dis(k) U(k)  (num_bands x num_wann), |w'_j> = sum_m V_mj |psi_m>,
  H^W(k) = V+(k) diag(eig(k)) V(k),  S_j^W(k) = V+(k) S_j^Bloch(k) V(k).
  The authored final gauge differs from the source bundle basis by a smooth
  rotation W(k): final operators (H_f, S_f, A_f) = W+ (H, S, A) W with the
  position connection law  A_f(k) = W+ A(k) W + i W+ d_k W  (d_k in 1/Ang;
  PLUS i under D_k H = dH/dk - i[A,H]).
  spn holds the spin in the ab-initio eigenbasis: S_E = U_eig+ S U_eig.
  r.dat is ALREADY final Wannier gauge: read directly, no second rotation;
  centers = diagonal on-site position entries, counted once.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)


def _kpts_grid(N):
    return np.array([[i / N, j / N, 0.0] for i in range(N) for j in range(N)])


def _fourier(kpts, R, mat):
    out = np.zeros((len(kpts), mat.shape[-2], mat.shape[-1]), complex)
    for ir, Rv in enumerate(R):
        ph = np.exp(2j * np.pi * (Rv[0] * kpts[:, 0] + Rv[1] * kpts[:, 1]))
        out += ph[:, None, None] * mat[ir]
    return out


def smooth_unitary(N, n, seed=7, rate=0.5):
    """Smooth k-dependent n x n unitary W(k) = expm(-i f(q) K0); nonidentity.
    Returns (W_grid, K0, rate) with W_grid reshaped (N*N, n, n) in grid_q
    order and f(q) = rate*(sin 2pi qx + cos 2pi qy). The analytic derivative
    is W+ dW/dq = -i (df/dq) K0 (K0 commutes with W), used for the exact
    position connection term (no finite-difference error)."""
    rng = np.random.default_rng(seed)
    K0 = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    K0 = K0 + K0.conj().T
    Nk = N * N
    W = np.zeros((Nk, n, n), complex)
    for k in range(Nk):
        qx, qy = (k // N) / N, (k % N) / N   # grid_q row-major (i outer, j inner)
        K = K0 * rate * (np.sin(2 * np.pi * qx) + np.cos(2 * np.pi * qy))
        ev, ve = np.linalg.eigh(K)
        W[k] = (ve * np.exp(-1j * ev)) @ ve.conj().T
    return W, K0, rate


def _window_R(w):
    return np.array([[i, j] for i in range(-w, w + 1) for j in range(-w, w + 1)],
                    dtype=int)


def author_w90_files(bundle, outdir, N, label="v2_positive",
                     gauge="identity", disentangle=False, seed=7,
                     rate=0.02, r_window=None):
    """Write official-syntax synthetic files reproducing `bundle` in the
    final Wannier gauge.

    gauge: 'identity'  -> final gauge == the bundle basis (W = I); the file
           V = U_eig (the exact diagonalizer) is still nonidentity.
           'smooth'    -> final gauge = smooth nonidentity W(k) with the
           exact analytic connection term (rate default 0.02 so the Fourier
           content of the k-dependent rotation fits the R window).
    disentangle: True -> 6 ab-initio bands, num_wann = 4,
           V = [U_eig+ W; 0] (6x4), U = W (4x4), U_dis = [U_eig+; 0] (6x4).
    r_window: half-width of the R window for the hr/r.dat R list; default
    None -> the bundle's own translations for the identity gauge, [-w, w]^2
    with w=5 for the smooth/disentangled gauges (captures the rotation's
    Fourier content).
    Returns (dir, receipt)."""
    D = Path(outdir)
    D.mkdir(parents=True, exist_ok=True)
    Nk = N * N
    kpts = _kpts_grid(N)
    if gauge == "identity":
        R = bundle.translations
    else:
        R = _window_R(r_window if r_window is not None else 5)
    Hq = _fourier(kpts, R if gauge == "identity" else bundle.translations,
                  bundle.h_r)
    Sq = np.stack([_fourier(kpts, bundle.translations, bundle.spin_r[:, j])
                   for j in range(3)], axis=1)                    # (Nk,3,4,4)
    Aq = np.stack([_fourier(kpts, bundle.translations, bundle.position_r[:, j])
                   for j in range(3)], axis=1)                    # (Nk,3,4,4)

    Ueig = np.array([np.linalg.eigh(Hq[k])[1] for k in range(Nk)])  # (Nk,4,4)
    K0 = None
    if gauge == "identity":
        Wf = np.tile(I4, (Nk, 1, 1))
    else:
        Wf, K0, _ = smooth_unitary(N, 4, seed=seed, rate=rate)
    a1, a2 = bundle.lattice[0, 0], bundle.lattice[1, 1]

    # final-gauge operators (joint rotation W with the position connection;
    # the connection term uses the EXACT analytic derivative W+ dW/dq =
    # -i (df/dq) K0, f = rate*(sin 2pi qx + cos 2pi qy))
    Hf = np.array([Wf[k].conj().T @ Hq[k] @ Wf[k] for k in range(Nk)])
    Sf = np.stack([np.stack([Wf[k].conj().T @ Sq[k, j] @ Wf[k]
                             for k in range(Nk)]) for j in range(3)], axis=1)
    Af = np.empty((Nk, 3, 4, 4), complex)
    for k in range(Nk):
        qx, qy = kpts[k, 0], kpts[k, 1]
        c1 = (a1 * rate * np.cos(2 * np.pi * qx) * K0 if K0 is not None else 0j)
        c2 = (-a2 * rate * np.sin(2 * np.pi * qy) * K0 if K0 is not None else 0j)
        for j in range(3):
            Af[k, j] = Wf[k].conj().T @ Aq[k, j] @ Wf[k]
        Af[k, 0] = Af[k, 0] + c1
        Af[k, 1] = Af[k, 1] + c2
    # Az: 2D slab -> no a3-derivative connection term (slab, Rz=0)

    # ab-initio eigenbasis spin (spn content)
    SE = np.stack([np.stack([Ueig[k].conj().T @ Sq[k, j] @ Ueig[k]
                             for k in range(Nk)]) for j in range(3)], axis=1)

    if disentangle:
        nb, nw = 6, 4
        eps = np.array([np.linalg.eigvalsh(Hq[k]) for k in range(Nk)])
        # V (6x4) = [U_eig+ W; 0];  U (4x4) = W;  U_dis (6x4) = V U+ = [U_eig+; 0]
        V_list = np.zeros((Nk, 6, 4), complex)
        V_list[:, :4, :] = np.stack([Ueig[k].conj().T @ Wf[k]
                                     for k in range(Nk)])
        U_list = Wf
        Ud_list = np.zeros((Nk, 6, 4), complex)
        Ud_list[:, :4, :] = np.stack([Ueig[k].conj().T for k in range(Nk)])
        # full 6x6 ab-initio spin: blockdiag(SE, sigma_far)
        SE6 = np.zeros((Nk, 3, 6, 6), complex)
        SE6[:, :, :4, :4] = SE
        for j in range(3):
            SE6[:, j, 4:, 4:] = np.array([SX, SY, SZ])[j]
        spn = SE6
        eig = np.concatenate([eps, np.full((Nk, 2), 25.0), np.full((Nk, 2), 26.0)],
                             axis=1)  # placeholder, overwritten below
        eig = np.column_stack([eps, np.full(Nk, 25.0), np.full(Nk, 26.0)])
    else:
        nb, nw = 4, 4
        eps = np.array([np.linalg.eigvalsh(Hq[k]) for k in range(Nk)])
        V_list = np.stack([Ueig[k].conj().T @ Wf[k] for k in range(Nk)])  # 4x4
        U_list = V_list          # u.mat: the full gauge rotation (V with U_dis=I)
        Ud_list = None
        spn = SE
        eig = eps

    # ---- hr.dat (final Wannier-gauge H(R); ndegen applied, deg all 1)
    H_R = _inverse_fourier(kpts, R, Hf)
    lines = [f" created on 23Oct2026 at 00:00:00 (synthetic fixture, {label})",
             f" {nw}", f" {len(R)}", " 1" * len(R)]
    for ir in range(len(R)):
        for m in range(nw):
            for n in range(nw):
                z = H_R[ir, m, n]
                lines.append(f"  {R[ir][0]:3d} {R[ir][1]:3d} 0 {m+1:2d} {n+1:2d}"
                             f" {z.real: .10f} {z.imag: .10f}")
    (D / "seed_hr.dat").write_text("\n".join(lines) + "\n")

    # ---- r.dat (final-gauge position connection; no degeneracy block)
    A_R = np.stack([_inverse_fourier(kpts, R, Af[:, j]) for j in range(3)],
                   axis=1)                                        # (nr,3,4,4)
    lines = [f" created on 23Oct2026 at 00:00:00 (synthetic fixture, {label})",
             f" {nw}", f" {len(R)}"]
    for ir in range(len(R)):
        for m in range(nw):
            for n in range(nw):
                row = [f"  {R[ir][0]:3d} {R[ir][1]:3d} 0 {m+1:2d} {n+1:2d}"]
                for j in range(3):
                    z = A_R[ir, j, m, n]
                    row.append(f" {z.real: .10f} {z.imag: .10f}")
                lines.append("".join(row))
    (D / "seed_r.dat").write_text("\n".join(lines) + "\n")

    # ---- spn (formatted, ab-initio eigenbasis spin, upper triangle;
    # loop order m=1..nb, n=1..m, i=x,y,z with i FASTEST per the docs)
    lines = [f" created on 23Oct2026 at 00:00:00 (synthetic spn, {label})",
             f" {nb} {Nk}"]
    for k in range(Nk):
        for m in range(1, nb + 1):
            for n in range(1, m + 1):
                for j in range(3):
                    z = spn[k, j][n - 1, m - 1]
                    lines.append(f" {z.real: .16f} {z.imag: .16f}")
    (D / "seed.spn").write_text("\n".join(lines) + "\n")

    # ---- u.mat / u_dis.mat / v.mat
    # Official Wannier90 U-family files store matrices in COLUMN-MAJOR order
    # (rows change fastest, then columns; per the files docs) - loop c outer,
    # r inner. A row-major file would silently transpose a genuine complex
    # nonsymmetric U (external-format golden test in test_import...).
    def mat_block(M_list, nrows, ncols):
        out = []
        for k in range(Nk):
            out.append(f"  {kpts[k,0]:.10f} {kpts[k,1]:.10f} {kpts[k,2]:.10f}")
            for c in range(ncols):
                for r in range(nrows):
                    z = M_list[k][r, c]
                    out.append(f" {z.real: .16f} {z.imag: .16f}")
            out.append("")
        return out

    uu = ([f" created on 23Oct2026 at 00:00:00 (synthetic u, {label})",
           f" {Nk} {nw} {nw}", ""] + mat_block(U_list, nw, nw)[:-1])
    (D / "seed_u.mat").write_text("\n".join(uu) + "\n")
    if disentangle:
        ud = ([f" created on 23Oct2026 at 00:00:00 (synthetic u_dis, {label})",
               f" {Nk} {nw} {nb}", ""] + mat_block(Ud_list, nb, nw)[:-1])
        (D / "seed_u_dis.mat").write_text("\n".join(ud) + "\n")
        vv = ([f" created on 23Oct2026 at 00:00:00 (synthetic v, {label})",
               f" {Nk} {nw} {nb}", ""] + mat_block(V_list, nb, nw)[:-1])
        (D / "seed_v.mat").write_text("\n".join(vv) + "\n")

    # ---- eig
    elines = [f" {b+1:3d} {k+1:5d} {eig[k, b]: .12f}"
              for k in range(Nk) for b in range(nb)]
    (D / "seed.eig").write_text("\n".join(elines) + "\n")

    # ---- win (explicit kpoints; orthorhombic 3x3; saxis; ndegen applied)
    a1v, a2v = bundle.lattice[0], bundle.lattice[1]
    a3v = np.array([0.0, 0.0, 20.0])
    win = (f"num_bands = {nb}\nnum_wann = {nw}\nspinors = .true.\n"
           f"saxis = 0 0 1\nwrite_ndegen_applied = .true.\n"
           f"use_ws_distance = .false.\nbegin unit_cell_cart\n"
           f" {a1v[0]:.10f} {a1v[1]:.10f} {a1v[2]:.10f}\n"
           f" {a2v[0]:.10f} {a2v[1]:.10f} {a2v[2]:.10f}\n"
           f" {a3v[0]:.10f} {a3v[1]:.10f} {a3v[2]:.10f}\n"
           f"end unit_cell_cart\nbegin kpoints\n")
    win += "".join(f"  {k[0]:.10f} {k[1]:.10f} {k[2]:.10f}\n" for k in kpts)
    win += "end kpoints\n"
    (D / "seed.win").write_text(win)

    receipt = {"label": label, "gauge": gauge, "disentangle": disentangle,
               "N": N, "Nk": Nk, "num_wann": nw, "num_bands": nb,
               "k_order": "grid_q row-major (i outer, j inner)",
               "final_gauge_rotation": "W(k): "
               + ("identity" if gauge == "identity" else "smooth nonidentity"),
               "files": sorted(p.name for p in D.iterdir())}
    return D, receipt


def _inverse_fourier(kpts, R, X_k):
    """X(R) = mean_q exp(-2pi i R.q) X(q) on the uniform grid."""
    out = []
    for Rv in R:
        ph = np.exp(-2j * np.pi * (Rv[0] * kpts[:, 0] + Rv[1] * kpts[:, 1]))
        out.append(np.mean(ph[:, None, None] * X_k, axis=0))
    return np.array(out)
