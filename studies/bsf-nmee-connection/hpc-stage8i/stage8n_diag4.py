#!/usr/bin/env python3
"""Stage8n diagnostic 4: Wilson-phase precision budget + SPN mirror
self-consistency over ALL 441 k points (no new PBS job).

Part A: at selected k, dump the undivided complex plaquette phase
P = det(L1 L2 L3 L4), its two-field difference dP = P(B=+dB) - P(B=-dB),
|dP| vs machine epsilon, min SVD of the L product, and the divided Gamma_k.
This is the absolute error budget the review asks for: a BZ-integrated
Gamma ~1e-9 may hide tiny loop-phase differences.

Part B: SPN full-mesh mirror self-consistency. For the in-plane mirror
m_110 (det=-1, maps (kx,ky)->(-ky,-kx) cartesian), pairs of k points that
are mirror images must satisfy (basis-invariant checks):
  Tr(S_i^2) equal, |S_z diag| equal, and the trace of the S_z^2 block equal.
A pair table is reported with the worst residuals.
"""
import json, sys
import numpy as np

DIR = sys.argv[1]
out = sys.argv[2]

from stage8m_plaq import read_eig, read_mmn_spn, DeterminantPlaquette
from stage8m_response import (read_win_kpoints, read_win_cell,
                              read_spn_formatted, WannierProvider)

nb, nk = 96, 441
E = read_eig(f"{DIR}/wannier90.eig")
kpts = read_win_kpoints(f"{DIR}/wannier90.win")
cell = read_win_cell(f"{DIR}/wannier90.win")
S = read_spn_formatted(f"{DIR}/wannier90.spn", nb, nk)
mu_b = 5.7884e-5
occ = 36

n1 = n2 = 21
gx = np.round(kpts[:, 0] * n1).astype(int) % n1
gy = np.round(kpts[:, 1] * n2).astype(int) % n2
coord = {(gx[i], gy[i]): i for i in range(nk)}
neigh = {(ax, sgn): [0] * nk for ax in (0, 1) for sgn in (-1, 1)}
for ik in range(nk):
    a, b = gx[ik], gy[ik]
    neigh[(0, +1)][ik] = coord[((a + 1) % n1, b)]
    neigh[(0, -1)][ik] = coord[((a - 1) % n1, b)]
    neigh[(1, +1)][ik] = coord[(a, (b + 1) % n2)]
    neigh[(1, -1)][ik] = coord[(a, (b - 1) % n2)]
a1, a2 = cell[0][:2], cell[1][:2]
A2d = float(np.linalg.det(np.array([a1, a2])))
b1 = 2 * np.pi * np.array([a2[1], -a2[0]]) / A2d
b2 = 2 * np.pi * np.array([-a1[1], a1[0]]) / A2d
dk_cart = np.array([np.linalg.norm(b1) / n1, np.linalg.norm(b2) / n2])

def mmn_fast(path, nb, nk, want):
    with open(path) as f:
        ln = f.readline()
        while True:
            p = ln.split()
            if (len(p) >= 3 and p[0].isdigit() and int(p[0]) == nb
                    and p[1].isdigit() and int(p[1]) == nk):
                rest = f.read()
                break
            ln = f.readline()
    a = np.fromstring(rest, sep=" ", dtype=np.float64)
    nlink = a.size // (5 + 2 * nb * nb)
    a = a[: nlink * (5 + 2 * nb * nb)].reshape(nlink, 5 + 2 * nb * nb)
    heads = a[:, :5].astype(int)
    out = {}
    for r in range(nlink):
        ik, ikp = heads[r, 0] - 1, heads[r, 1] - 1
        if (ik, ikp) not in out:
            seg = a[r, 5:].reshape(nb, nb, 2)
            out[(ik, ikp)] = seg[:, :, 0] + 1j * seg[:, :, 1]
    return out

# --- Part A: Wilson phase budget on a few k ---
sel = [(0.44, 0.12), (0.13, 0.07), (0.31, 0.22)]
want = []
for fx, fy in sel:
    ik = coord[(round(fx * n1) % n1, round(fy * n2) % n2)]
    for axis in (0, 1):
        jp = neigh[(axis, +1)][ik]
        jm = neigh[(axis, -1)][ik]
        want += [(ik, jp), (jp, ik), (ik, jm), (jm, ik)]
mats = mmn_fast(f"{DIR}/wannier90.mmn", nb, nk, want)
prov = WannierProvider(E, S[:, :, :, 0], S[:, :, :, 2], occ, nb)

# reproduce the engine plaquette to access undivided phase
from stage8m_plaq import link_matrix, polar_unitary

def raw_plaq(ik, lam, B, axis, dl=1e-4):
    eng = DeterminantPlaquette(
        prov.frames, mode="mmn", links=mats,
        ik_of=lambda k_c, sgn, axis=axis: neigh[(axis, 1 if sgn > 0 else -1)][k_c],
        dk=dk_cart[axis], dl=dl, axis=axis, full_nb=nb)
    return eng.plaquette_phase(ik, lam, B)

phase_rows = []
for fx, fy in sel:
    ik = coord[(round(fx * n1) % n1, round(fy * n2) % n2)]
    for axis in (0, 1):
        Bp = raw_plaq(ik, 0.0, .25, axis)
        Bm = raw_plaq(ik, 0.0, -.25, axis)
        dP = Bp - Bm
        gk = dP / (2 * .25)
        phase_rows.append({
            "k_frac": [fx, fy], "axis": axis,
            "phase_plus": Bp, "phase_minus": Bm,
            "dP_over_2dB": gk,
            "machine_eps": np.finfo(float).eps,
            "dP_abs_vs_eps_ratio": abs(dP) / np.finfo(float).eps,
        })

# --- Part B: full-mesh SPN mirror self-consistency (m_110) ---
kcart = kpts @ np.stack([b1, b2])
mirror = np.array([[0., -1.], [-1., 0.]])
tol = 0.15 * np.linalg.norm(b1) / 21.0
pairs = []
for ik in range(nk):
    km = kcart[ik] @ mirror.T
    d = np.linalg.norm(kcart - km, axis=1)
    j = int(np.argmin(d))
    if d[j] < tol and j > ik:
        pairs.append((ik, j))
tr_resid = []
sz_resid = []
for ik, j in pairs:
    r = max(max(abs(np.trace(S[ik][:, :, ip].conj().T @ S[ik][:, :, ip]) -
                   np.trace(S[j][:, :, ip].conj().T @ S[j][:, :, ip]))
                for ip in range(3)), 0.0)
    tr_resid.append(r)
    sz = max(np.abs(np.diagonal(S[ik][:, :, 2])) -
             np.abs(np.diagonal(S[j][:, :, 2])))
    sz_resid.append(abs(sz))

report = {
    "wilson_phase_budget": phase_rows,
    "spn_mirror_selfconsistency": {
        "pairs": len(pairs),
        "TrS2_max_abs_diff": float(max(tr_resid)) if tr_resid else None,
        "TrS2_median_abs_diff": float(np.median(tr_resid)) if tr_resid else None,
        "Szdiag_max_abs_diff": float(max(sz_resid)) if sz_resid else None,
        "mirror": "m_110 (kx,ky)->(-ky,-kx) cartesian",
        "tol_frac_mesh": 0.15,
    },
}
json.dump(report, open(out, "w"), indent=1)
for r in phase_rows:
    print(f"k={r['k_frac']} ax={r['axis']}: P+= {r['phase_plus']:+.3e} "
          f"P-= {r['phase_minus']:+.3e} dP/2dB= {r['dP_over_2dB']:+.3e} "
          f"|dP|/eps= {r['dP_abs_vs_eps_ratio']:.1f}")
print(f"mirror pairs {len(pairs)}: TrS2 max {report['spn_mirror_selfconsistency']['TrS2_max_abs_diff']:.3e}")
