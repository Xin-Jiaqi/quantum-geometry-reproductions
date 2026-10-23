#!/usr/bin/env python3
"""Stage8n diagnostic 3: independent Kubo/projector B-derivative cross-check
on the real 96-band outputs (no new PBS job).

Per selected k points:
  F_kubo(B) = sum_{n occ, m unocc} 2 Im[ <n|S_i/2|m> <m|dH/dk_j|n> ] / (E_n-E_m)^2
  Gamma_k = [F_kubo(+dB) - F_kubo(-dB)] / (2 dB)
with H(B) = diag E - muB*B*S_z (frozen), band velocities from central overlap
gradients of the committed MMN neighbour links, SPN/2 spin source, full 36
occupied bands. Compared with the DeterminantPlaquette engine Gamma_k at the
same k (B=0 central diff, dB=0.25).
"""
import json, sys
import numpy as np

DIR = sys.argv[1]
out = sys.argv[2]

from stage8m_plaq import read_eig, read_mmn_spn, link_matrix, DeterminantPlaquette
from stage8m_response import read_win_kpoints, read_win_cell, read_spn_formatted, WannierProvider

nb, nk = 96, 441
E = read_eig(f"{DIR}/wannier90.eig")
kpts = read_win_kpoints(f"{DIR}/wannier90.win")
cell = read_win_cell(f"{DIR}/wannier90.win")
S = read_spn_formatted(f"{DIR}/wannier90.spn", nb, nk)
def read_mmn_fast(path, nb, nk, want):
    """Vectorized reader: returns dict {(ik,ikp): M} only for requested
    (ik,ikp) bonds. Block layout per link: 5 header tokens then nb*nb pairs.
    Skips the VASP comment line(s) before the numeric header."""
    with open(path) as f:
        ln = f.readline()
        while True:
            p = ln.split()
            if (len(p) >= 3 and p[0].isdigit() and int(p[0]) == nb
                    and p[1].isdigit() and int(p[1]) == nk):
                rest = f.read()   # data only; header line consumed
                break
            ln = f.readline()
            if not ln:
                raise ValueError("no numeric header in mmn")
    a = np.fromstring(rest, sep=" ", dtype=np.float64)
    nlink = a.size // (5 + 2 * nb * nb)
    a = a[: nlink * (5 + 2 * nb * nb)].reshape(nlink, 5 + 2 * nb * nb)
    heads = a[:, :5].astype(int)
    want_set = set(want)
    out = {}
    for r in range(nlink):
        ik, ikp = heads[r, 0] - 1, heads[r, 1] - 1
        if (ik, ikp) in want_set and (ik, ikp) not in out:
            seg = a[r, 5:].reshape(nb, nb, 2)
            out[(ik, ikp)] = seg[:, :, 0] + 1j * seg[:, :, 1]
    # missing bonds are fine: callers fall back to the conjugate transpose
    # of the reverse bond (link_matrix convention)
    return out

mu_b = 5.7884e-5
occ = 36
# selected k: several low-symmetry points away from Gamma (all have the same
# 1.153 eV occ-unocc gap; avoid exact high-symmetry for independent check)
sel_frac = [(0.13, 0.07), (0.31, 0.22), (0.05, 0.41), (0.44, 0.12), (0.24, 0.33)]

# k-graph neighbours (periodic index lookup, same as response pipeline)
n1 = n2 = 21
gx = np.round(kpts[:, 0] * n1).astype(int) % n1
gy = np.round(kpts[:, 1] * n2).astype(int) % n2
coord = {(gx[i], gy[i]): i for i in range(nk)}
if len(coord) != nk:
    raise ValueError("win kpoints do not cover the mesh")
neigh = {(ax, sgn): [0] * nk for ax in (0, 1) for sgn in (-1, 1)}
for ik in range(nk):
    a, b = gx[ik], gy[ik]
    neigh[(0, +1)][ik] = coord[((a + 1) % n1, b)]
    neigh[(0, -1)][ik] = coord[((a - 1) % n1, b)]
    neigh[(1, +1)][ik] = coord[(a, (b + 1) % n2)]
    neigh[(1, -1)][ik] = coord[(a, (b - 1) % n2)]
want = []
for fx, fy in sel_frac:
    ik = coord[(round(fx * n1) % n1, round(fy * n2) % n2)]
    for axis in (0, 1):
        jp = neigh[(axis, +1)][ik]
        jm = neigh[(axis, -1)][ik]
        want += [(ik, jp), (jp, ik), (ik, jm), (jm, ik)]
mats = read_mmn_fast(f"{DIR}/wannier90.mmn", nb, nk, want)
a1, a2 = cell[0][:2], cell[1][:2]
A2d = float(np.linalg.det(np.array([a1, a2])))
b1 = 2 * np.pi * np.array([a2[1], -a2[0]]) / A2d
b2 = 2 * np.pi * np.array([-a1[1], a1[0]]) / A2d
dk_cart = np.array([np.linalg.norm(b1) / n1, np.linalg.norm(b2) / n2])

def F_kubo(ik, axis, B, src=0):
    """F(B) via independent Kubo/projector resolvent (36 occ bands).

    Gauge-consistent construction: all operators live in the VASP band basis
    (E diagonal, SPN/MMN as written); the Zeeman-perturbed eigenbasis at B is
    u(B) = U(B) (eigh of H(B) in the VASP basis), and both the spin source and
    the velocity gradient are rotated by the SAME U(B), so the off-diagonal
    Im[<m|S_i|n><n|v|m>] is gauge invariant."""
    H = np.diag(E[ik]).astype(complex) - mu_b * B * S[ik][:, :, 2]
    e, U = np.linalg.eigh(H)
    # central overlap gradient along axis using real MMN links (VASP basis)
    jp = neigh[(axis, +1)][ik]
    jm = neigh[(axis, -1)][ik]
    Mp = mats[(ik, jp)] if (ik, jp) in mats else mats[(jp, ik)].conj().T
    Mm = mats[(ik, jm)] if (ik, jm) in mats else mats[(jm, ik)].conj().T
    dk = dk_cart[axis]
    dM = (Mp - Mm) / (2 * dk)              # <u_m(k)|d u_n(k)/dk> in VASP basis
    dM_b = U.T.conj() @ dM @ U            # rotate to u(B) basis
    S2_b = U.T.conj() @ (S[ik][:, :, src] / 2.0) @ U
    gam = 0.0
    for n in range(occ):
        for m in range(occ, nb):
            vnm = (e[n] - e[m]) * dM_b[m, n]
            d = (e[n] - e[m]) ** 2
            if d > 1e-12:
                gam += 2.0 * np.imag(S2_b[m, n] * vnm) / d
    return float(gam)

rows = []

for fx, fy in sel_frac:
    ik = coord[(round(fx * n1) % n1, round(fy * n2) % n2)]
    g0 = float(np.min(E[ik][occ] - E[ik][occ - 1]))
    for axis in (0, 1):
        prov = WannierProvider(E, S[:, :, :, 0], S[:, :, :, 2], occ, nb)
        eng = DeterminantPlaquette(
            prov.frames, mode="mmn", links=mats,
            ik_of=lambda k_c, sgn, axis=axis: neigh[(axis, 1 if sgn > 0 else -1)][k_c],
            dk=dk_cart[axis], dl=1e-4, axis=axis, full_nb=nb)
        g_plaq = eng.Gamma_k(ik, 0.0, 0.0, .25)
        g_kub = (F_kubo(ik, axis, .25) - F_kubo(ik, axis, -.25)) / 0.5
        rows.append({"k_frac": [fx, fy], "axis": axis, "gap_eV": g0,
                     "plaq": g_plaq, "kubo": g_kub,
                     "rel": abs(g_kub - g_plaq) / max(abs(g_plaq), 1e-30)})
json.dump({"rows": rows}, open(out, "w"), indent=1)
for r in rows:
    print(f"k={r['k_frac']} ax={r['axis']}: plaq={r['plaq']:+.4e} kubo={r['kubo']:+.4e} rel={r['rel']:.2e}")
