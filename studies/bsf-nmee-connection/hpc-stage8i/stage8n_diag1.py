#!/usr/bin/env python3
"""Stage8n diagnostic 1: window closure + SPN symmetry self-consistency.

Runs on the Japan IMR login node against the committed 96-band
wannier90.{eig,spn,win} outputs (no new PBS job, no raw upload).

Outputs (JSON):
 - band-96 window profile: min/max E96 over k, min gap E96-E95, k-location
 - SPN self-consistency: per-k S_z diag bounds, Tr(S_i^2) invariance under
   in-plane mirror-equivalent k pairs (mirror ops of the stack), and the
   k->k' pair list used.
 - window-closure note: with NBANDS=96 no band-97 energies exist in the
   committed outputs, so band-96/97 crossings cannot be checked from data;
   we instead bound the closure risk via max E96 and the min E96-E95 gap.
"""
import json, sys
from pathlib import Path
import numpy as np

DIR = sys.argv[1]
out = sys.argv[2]

def read_eig(path):
    nk, nb = 441, 96
    E = np.zeros((nk, nb))
    for ln in open(path):
        p = ln.split()
        if len(p) < 3:
            continue
        b = int(p[0]) - 1
        k = int(p[1]) - 1
        E[k, b] = float(p[2])
    return E

def read_win_kpoints(path):
    kpts = []
    lines = open(path).read().splitlines()
    i = lines.index("begin kpoints")
    for ln in lines[i + 1:]:
        if ln.strip().lower().startswith("end"):
            break
        p = ln.split()
        if len(p) >= 2:
            kpts.append([float(p[0]), float(p[1])])
    return np.array(kpts)

def read_spn_sample(path, nk, nb, stride=17):
    """Read SPN (nb nk header, 3*Ntri rows per k, Re Im)."""
    Ntri = nb * (nb + 1) // 2
    rows = 3 * Ntri
    indn, indm = np.tril_indices(nb)
    S = {}
    with open(path) as f:
        f.readline()  # 'Written by ...' comment line
        hdr = f.readline().split()
        nnb, nnk = map(int, hdr[:2])
        assert nnb == nb and nnk == nk, (nnb, nnk)
        for ik in range(nk):
            A3 = np.zeros((3, nb, nb), complex)
            for r in range(rows):
                z = complex(*map(float, f.readline().split()))
                tri = r // 3
                ip = r % 3
                A3[ip, indn[tri], indm[tri]] = z
            A = np.zeros((nb, nb, 3), complex)
            for ip in range(3):
                A[:, :, ip] = A3[ip] + np.tril(A3[ip], -1).conj().T
            if ik % stride == 0:
                S[ik] = A
    return S

nb, nk = 96, 441
E = read_eig(f"{DIR}/wannier90.eig")
kpts = read_win_kpoints(f"{DIR}/wannier90.win")
S = read_spn_sample(f"{DIR}/wannier90.spn", nk, nb)

# --- window closure profile ---
E96 = E[:, 95]
E95 = E[:, 94]
gap9596 = E96 - E95
k96max = int(np.argmax(E96))
k96min = int(np.argmin(E96))

# --- SPN self-consistency under in-plane mirror-equivalent k pairs ---
# mirror ops of the cell (from spglib analysis): GLIDE: R=[[0,-1],[-1,0]]
# (m_110, det=-1) maps k=(kx,ky)->(-ky,-kx) in cartesian; AB: m_110 and m_010.
# In fractional (a*,b*) coords the mirror acts as the transpose of the
# cartesian rotation conjugated by the reciprocal lattice; simplest robust
# check: compare Tr(S_i^2) (basis-invariant) over k points related by the
# in-plane mirror, using cartesian k.
import itertools
a1 = np.array([3.3190066244679999, 0.0])
a2 = np.array([-1.6595033122339999, 2.8743440521179999])
A2d = a1[0]*a2[1] - a1[1]*a2[0]
b1 = 2*np.pi*np.array([a2[1], -a2[0]])/A2d
b2 = 2*np.pi*np.array([-a1[1], a1[0]])/A2d
kcart = kpts @ np.stack([b1, b2])
mirrors = [np.array([[0., -1.], [-1., 0.]]),   # m_110
           np.array([[1., 0.], [0., -1.]])]    # m_010 (if present)
tol = 1.0/21 * 0.15  # fractional tolerance of half a mesh step
def find_k(kc, exclude):
    for i in range(len(kpts)):
        if i in exclude:
            continue
        if np.linalg.norm(kcart[i] - kc) < tol * np.linalg.norm(b1):
            return i
    return None
spn_tr = {}
tr_resid = {}
for ik, A in S.items():
    spn_tr[ik] = [float(np.trace(A[:, :, ip].conj().T @ A[:, :, ip]).real)
                  for ip in range(3)]
for mk in mirrors:
    n = 0
    worst = 0.0
    for ik, A in S.items():
        km = kcart[ik] @ mk.T
        j = find_k(km, {ik})
        if j is not None and j in S:
            d = max(abs(spn_tr[ik][ip] - spn_tr[j][ip]) for ip in range(3))
            worst = max(worst, d)
            n += 1
    tr_resid[str(mk.tolist())] = {"pairs": n, "max_TrDiff": float(worst)}

# --- Hermiticity & Pauli bounds on sampled k ---
herm_max = 0.0
sz_diag_max = 0.0
for ik, A in S.items():
    herm_max = max(herm_max, np.max(np.abs(A.transpose(1, 0, 2).conj() - A)))
    sz_diag_max = max(sz_diag_max, np.max(np.abs(np.diagonal(A[:, :, 2]))))

report = {
    "window_closure": {
        "E96_min_eV": float(E96.min()),
        "E96_max_eV": float(E96.max()),
        "E96_at_max_k": int(k96max),
        "E96_at_min_k": int(k96min),
        "min_gap_E95_E96_eV": float(gap9596.min()),
        "k_of_min_gap": int(np.argmin(gap9596)),
        "note": ("NBANDS=96: no band-97 energies in committed wannier90.eig, "
                 "so 96/97 crossings cannot be verified from data; closure "
                 "risk bounded by max E96 and min E95-E96 gap."),
    },
    "spn_sampling": {
        "kpoints_sampled": len(S),
        "stride": 17,
        "hermiticity_max_dev": float(herm_max),
        "S_z_diag_max": float(sz_diag_max),
    },
    "spn_tr_invariance": {
        "per_k_TrS2": {str(ik): spn_tr[ik] for ik in sorted(S)},
        "mirror_residuals": tr_resid,
    },
}
json.dump(report, open(out, "w"), indent=1)
print(f"wrote {out}: E96 [{E96.min():.3f},{E96.max():.3f}] eV, min E95-E96 gap {gap9596.min():.4f} eV")
