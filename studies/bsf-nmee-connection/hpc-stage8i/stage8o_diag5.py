#!/usr/bin/env python3
"""Stage8o real-data Kubo three-way cross-check (zero new PBS jobs).

At selected well-gapped low-symmetry k of the 96-band AB stack:
  (1) resolvent Kubo:  F = 2 sum Im[S2_zb[n,m] dM_cov[m,n]]/(E_n-E_m)
      with gauge-covariant dM_cov = (U(k)^+ M(k,k+) U(k+) - U(k)^+ M(k,k-) U(k-))/(2 dk)
  (2) projector finite-difference: F = i Tr P[dP/dl, dP/dk] (P in u(B) basis)
  (3) Wilson plaquette engine Gamma_k (B=0 central diff, dB=0.25)
Also reports the raw loop angle phi = |arg det(L1L2L3L4)| and its two-field
difference for the Wilson-phase precision budget.
"""
import json, sys
import numpy as np

DIR = sys.argv[1]
out = sys.argv[2]

from stage8m_plaq import read_eig, read_mmn_spn, DeterminantPlaquette
from stage8m_response import read_win_kpoints, read_win_cell, read_spn_formatted, WannierProvider
from stage8o_kubo_cross import resolvent_gamma_k, projector_gamma_k, dM_covariant

nb, nk = 96, 441
mu_b = 5.7884e-5
occ = 36
E = read_eig(f"{DIR}/wannier90.eig")
kpts = read_win_kpoints(f"{DIR}/wannier90.win")
cell = read_win_cell(f"{DIR}/wannier90.win")
S = read_spn_formatted(f"{DIR}/wannier90.spn", nb, nk)

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

def bond(ik, j):
    if (ik, j) in mats:
        return mats[(ik, j)]
    if (j, ik) in mats:
        return mats[(j, ik)].conj().T
    return None

# selected exact mesh k points (21x21 grid indices, reported as grid coords)
sel = [(8, 2), (7, 9)]
want = []
for i, j in sel:
    ik = coord[(i % n1, j % n2)]
    for axis in (0, 1):
        jp = neigh[(axis, +1)][ik]
        jm = neigh[(axis, -1)][ik]
        want += [(ik, jp), (jp, ik), (ik, jm), (jm, ik)]
mats = mmn_fast(f"{DIR}/wannier90.mmn", nb, nk, want)
prov = WannierProvider(E, S[:, :, :, 0], S[:, :, :, 2], occ, nb)

def H_B(ik, B):
    return np.diag(E[ik]).astype(complex) - mu_b * B * S[ik][:, :, 2]

rows = []
try:
  for i, j in sel:
    ik = coord[(i % n1, j % n2)]
    g0 = float(np.min(E[ik][occ] - E[ik][occ - 1]))
    kfrac = [i / n1, j / n2]
    for axis in (0, 1):
        jp = neigh[(axis, +1)][ik]
        jm = neigh[(axis, -1)][ik]
        dk = dk_cart[axis]
        # --- plaquette engine ---
        eng = DeterminantPlaquette(
            prov.frames, mode="mmn", links=mats,
            ik_of=lambda k_c, sgn, axis=axis: neigh[(axis, 1 if sgn > 0 else -1)][k_c],
            dk=dk, dl=1e-4, axis=axis, full_nb=nb)
        g_plaq = eng.Gamma_k(ik, 0.0, 0.0, .25)
        # raw loop angle (undivided) at B=+dB
        phi_plus = abs(eng.plaquette_phase(ik, 0.0, .25)) * dk * 1e-4
        phi_minus = abs(eng.plaquette_phase(ik, 0.0, -.25)) * dk * 1e-4
        # --- resolvent + projector at B=+dB and -dB ---
        def kubo_at(B):
            eB, U = np.linalg.eigh(H_B(ik, B))
            Mp = bond(ik, jp)
            Mm = bond(ik, jm)
            if Mp is None or Mm is None:
                raise RuntimeError(f"missing MMN bond near k={ik} ax={axis}")
            _, Up = np.linalg.eigh(H_B(jp, B))
            _, Um = np.linalg.eigh(H_B(jm, B))
            dM = (U.T.conj() @ Mp @ Up - U.T.conj() @ Mm @ Um) / (2 * dk)
            S2_zb = U.T.conj() @ (S[ik][:, :, 0] / 2.0) @ U
            g_res, gmin = resolvent_gamma_k(eB, S2_zb, dM, occ)
            g_proj, _ = projector_gamma_k(E[ik], S[ik][:, :, 2], S[ik][:, :, 0],
                                          Mp, Mm, occ, dk, 1e-4, B)
            return g_res, g_proj
        gp, pp = kubo_at(.25)
        gm, pm = kubo_at(-.25)
        gk_res = (gp - gm) / 0.5
        gk_proj = (pp - pm) / 0.5
        rows.append({
            "k_grid": kfrac, "k_index": ik, "axis": axis,
            "min_occ_unocc_gap_eV": g0,
            "plaq_Gamma_k": g_plaq,
            "resolvent_Gamma_k": gk_res,
            "projector_Gamma_k": gk_proj,
            "rel_res_plaq": abs(gk_res - g_plaq) / max(abs(g_plaq), 1e-30),
            "rel_proj_res": abs(gk_proj - gk_res) / max(abs(gk_res), 1e-30),
            "raw_loop_angle_plus_rad": phi_plus,
            "raw_loop_angle_minus_rad": phi_minus,
            "raw_angle_diff_over_eps": abs(phi_plus - phi_minus) / np.finfo(float).eps,
        })
except Exception:
    import traceback
    traceback.print_exc()
    raise
json.dump({"rows": rows}, open(out, "w"), indent=1)
for r in rows:
    print(f"k={r['k_grid']} ax={r['axis']}: plaq={r['plaq_Gamma_k']:+.4e} "
          f"res={r['resolvent_Gamma_k']:+.4e} proj={r['projector_Gamma_k']:+.4e} "
          f"relRP={r['rel_res_plaq']:.2e} relPR={r['rel_proj_res']:.2e} "
          f"phi={r['raw_loop_angle_plus_rad']:.3e} dphi/eps={r['raw_angle_diff_over_eps']:.1e}")
