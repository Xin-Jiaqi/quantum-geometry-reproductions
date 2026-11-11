#!/usr/bin/env python3
"""Stage8m: coarse frozen-KS static spin-Zeeman response estimate from VASP
vasp2wannier90 postprocess files (wannier90.mmn/.spn/.eig/.win).

Physics (per STAGE8M taskbook): Gamma^{stat}_{iz;j} = d/dB_z of the mixed
parameter-space Berry curvature F_{lam_i,k_j} of the occupied projector,
with H(k,lam,B) = diag(E(k)) + lam*S_i/2 - muB*B*S_z (g=2, spin-only
Zeeman, frozen KS bands; S from wannier90.spn which stores Pauli-matrix
elements <u_m|sigma_ipol|u_n> -- verified |diag|>0.5 so it is sigma, not
sigma/2). NOT self-consistent magneto/electric response, NOT driven-DC
NMEE, NOT Edelstein.

Usage: stage8m_response.py --dir DIR [--windows 36,48,64,80,96]
       [--Blist 0.125,0.25] [--dllist 1e-4,2e-4,4e-4] [--source x,y]
       [--out out.json] [--kmesh 21,21] [--occ 36]
k points come from the .win kpoints block (fractional); mesh neighbours are
found by periodic index lookup. dk along a Cartesian axis is |b_axis|/N.
Full-BZ estimator: Gamma = (A_BZ/(2pi)^2) * (1/N_k) * sum_k Gamma(k).
"""
import argparse, json, time, os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np

from stage8m_plaq import (read_eig, read_mmn_spn, DeterminantPlaquette)

MU_B = 5.7883818060e-5   # eV/T


def parse_block(txt, begin, end):
    kl = txt.lower()
    i = kl.find(begin)
    j = kl.find(end)
    if i < 0 or j < 0:
        raise ValueError(f"block {begin!r}..{end!r} not found in .win")
    i += len(begin)
    return txt[i:j]


def read_win_kpoints(path):
    txt = Path(path).read_text()
    rows = []
    for line in parse_block(txt, "begin kpoints", "end kpoints").splitlines():
        p = [float(x) for x in line.split()]
        if len(p) >= 2:
            rows.append(p[:2])
    return np.array(rows)


def read_win_cell(path):
    txt = Path(path).read_text()
    vecs = []
    for line in parse_block(txt, "begin unit_cell_cart",
                            "end unit_cell_cart").splitlines():
        p = [float(x) for x in line.split()]
        if len(p) == 3:
            vecs.append(p)
    return np.array(vecs)          # rows are lattice vectors (Ang)


def read_spn_formatted(path, nb, nk):
    """wannier90.spn formatted: text header, 'nb nk', then per k-point
    3*(nb*(nb+1)//2) lines 'Re Im' (Fortran-ordered triangle per component).
    Returns data[ik, m, n, ipol] = <u_m|sigma_ipol|u_n> (Hermitian)."""
    Ntri = nb * (nb + 1) // 2
    indn, indm = np.tril_indices(nb)
    data = np.zeros((nk, nb, nb, 3), complex)
    with open(path) as f:
        f.readline()
        hdr = f.readline().split()
        if not (len(hdr) >= 2 and int(hdr[0]) == nb and int(hdr[1]) == nk):
            raise ValueError(f"spn header mismatch: {hdr}")
        for ik in range(nk):
            tmp = np.empty(3 * Ntri, complex)
            for t in range(3 * Ntri):
                p = f.readline().split()
                if len(p) < 2:
                    raise ValueError(f"short spn row at ik={ik} t={t}")
                tmp[t] = complex(float(p[0]), float(p[1]))
            A = np.zeros((3, nb, nb), complex)
            A[:, indn, indm] = tmp.reshape(3, Ntri, order="F")
            A[:, indm, indn] = A[:, indn, indm].conj()
            data[ik] = A.transpose(1, 2, 0)
    return data


def build_mesh(n1, n2):
    """Unused legacy helper kept for compatibility; neighbour tables are now
    built from the win kpoint coordinates (see run_response)."""
    return None, None


def build_hamiltonian(Ek, Sx, Sz, lam, B, window):
    H = np.diag(Ek[:window]).astype(complex)
    H[:window, :window] += lam * 0.5 * Sx[:window, :window]
    H[:window, :window] -= MU_B * B * Sz[:window, :window]
    return H


class WannierProvider:
    def __init__(self, E, Sx, Sz, occ, window, gap_min=0.0):
        self.E, self.Sx, self.Sz = E, Sx, Sz
        self.occ, self.window = occ, window
        self.gap_min = gap_min
        self.last_gap = None

    def frames(self, ik, lam, B):
        H = build_hamiltonian(self.E[ik], self.Sx[ik], self.Sz[ik],
                              lam, B, self.window)
        e, v = np.linalg.eigh(H)
        if self.occ < self.window:
            gap = float(e[self.occ] - e[self.occ - 1])
            if gap <= self.gap_min:
                raise ValueError(f"gap closed at ik={ik} lam={lam} B={B}: {gap}")
            self.last_gap = gap
        return np.ascontiguousarray(v[:, :self.occ])


_GLOBAL = {}


def _compute_one(args_shared, source, window, B, dl, axis):
    """Compute Gamma_BZ for one (source, window, B, dl, axis) config."""
    E, S, mats, neigh, dk_cart, weight, occ, nk0, nb0, gap_min = args_shared
    Sx = S[:, :, :, 0 if source == "x" else 1]
    Sz = S[:, :, :, 2]
    prov = WannierProvider(E, Sx, Sz, occ, window, gap_min)
    eng = DeterminantPlaquette(
        prov.frames, mode="mmn", links=mats,
        ik_of=lambda k_c, sgn, axis=axis: neigh[(axis, 1 if sgn > 0 else -1)][k_c],
        dk=dk_cart[axis], dl=dl, axis=axis, full_nb=nb0)
    gsum = 0.0
    ming = float("inf")
    for ik in range(nk0):
        gsum += eng.Gamma_k(ik, 0.0, 0.0, B)
        if prov.last_gap is not None:
            ming = min(ming, prov.last_gap)
    mins = min(eng.min_svals) if eng.min_svals else float("inf")
    return {
        "Gamma_mesh_mean": float(gsum / nk0),
        "Gamma_BZ_Ang2": float(weight * gsum),
        "min_gap_eV": float(ming),
        "min_sval": float(mins),
    }


def run_response(dirp, *, windows, blist, dllist, sources, occ, kmesh,
                 jobs=1):
    dirp = Path(dirp)
    t0 = time.time()
    E = read_eig(dirp / "wannier90.eig")
    nb0, nk0, nn0, mats = read_mmn_spn(dirp / "wannier90.mmn")
    kfrac = read_win_kpoints(dirp / "wannier90.win")
    cell = read_win_cell(dirp / "wannier90.win")
    if kfrac.shape[0] != nk0:
        raise ValueError(f"win kpoints {kfrac.shape[0]} != mmn nk {nk0}")
    if nb0 != E.shape[1]:
        raise ValueError("mmn nb != eig nb")
    n1, n2 = kmesh
    if n1 * n2 != nk0:
        raise ValueError(f"kmesh {kmesh} != nk {nk0}")
    # fractional -> grid coords (VASP writes symmetric-order kpoints, so we
    # must map by coordinates, not assume any linear mesh ordering)
    gx = np.round(kfrac[:, 0] * n1).astype(int) % n1
    gy = np.round(kfrac[:, 1] * n2).astype(int) % n2
    coord = {}
    for ik in range(nk0):
        c = (gx[ik], gy[ik])
        if c in coord:
            raise ValueError(f"duplicate grid coord {c} in win kpoints")
        coord[c] = ik
    if len(coord) != nk0:
        raise ValueError("win kpoints do not cover the mesh")
    neigh = {(ax, sgn): [0] * nk0 for ax in (0, 1) for sgn in (-1, 1)}
    for ik in range(nk0):
        a, b = gx[ik], gy[ik]
        neigh[(0, +1)][ik] = coord[((a + 1) % n1, b)]
        neigh[(0, -1)][ik] = coord[((a - 1) % n1, b)]
        neigh[(1, +1)][ik] = coord[(a, (b + 1) % n2)]
        neigh[(1, -1)][ik] = coord[(a, (b - 1) % n2)]

    a1, a2 = cell[0][:2], cell[1][:2]
    A2d = float(np.linalg.det(np.array([a1, a2])))            # Ang^2
    b1 = 2 * np.pi * np.array([a2[1], -a2[0]]) / A2d
    b2 = 2 * np.pi * np.array([-a1[1], a1[0]]) / A2d
    ABZ = abs(b1[0] * b2[1] - b1[1] * b2[0])                  # Ang^-2
    dk_cart = np.array([np.linalg.norm(b1) / n1,
                        np.linalg.norm(b2) / n2])             # Ang^-1
    weight = ABZ / ((2 * np.pi) ** 2) / nk0

    t1 = time.time()
    S = read_spn_formatted(dirp / "wannier90.spn", nb0, nk0)
    t2 = time.time()
    herm = max(np.max(np.abs(S[ik].transpose(1, 0, 2).conj() - S[ik]))
               for ik in range(0, nk0, 17))
    if herm > 1e-4:
        raise ValueError(f"SPN not Hermitian: {herm}")
    szd = np.max(np.abs(np.diagonal(S, axis1=1, axis2=2)[:, :, 2]))
    if szd > 1.0 + 1e-6:
        raise ValueError(f"S_z diag exceeds Pauli bound: {szd}")

    results = {}
    tasks = []
    for source in sources:
        for window in windows:
            for B in blist:
                for dl in dllist:
                    tasks.append((source, window, B, dl))
    shared = (E, S, mats, neigh, dk_cart, weight, occ, nk0, nb0, 0.0)
    t3 = time.time()
    if jobs <= 1:
        for (source, window, B, dl) in tasks:
            per_axis = {}
            for axis in (0, 1):
                per_axis[axis] = _compute_one(shared, source, window, B, dl,
                                              axis)
            key = f"src{source}_win{window}_B{B}_dl{dl}"
            results[key] = {
                "axis_x": per_axis[0],
                "axis_y": per_axis[1],
                "sym_residual_abs": abs(per_axis[0]["Gamma_BZ_Ang2"]
                                        + per_axis[1]["Gamma_BZ_Ang2"]),
            }
            print(key, results[key]["axis_x"]["Gamma_BZ_Ang2"],
                  results[key]["axis_y"]["Gamma_BZ_Ang2"],
                  "gap", per_axis[0]["min_gap_eV"], flush=True)
    else:
        # one task per (config, axis) for good load balancing; threads share
        # the read-only arrays and numpy releases the GIL in LAPACK calls,
        # so a thread pool avoids fork/OpenBLAS deadlocks seen with process
        # pools on this cluster
        axis_tasks = [(s, w, B, dl, ax)
                      for (s, w, B, dl) in tasks for ax in (0, 1)]
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futures = {ex.submit(_compute_one, shared, *t): t
                       for t in axis_tasks}
            for fut in futures:
                t = futures[fut]
                source, window, B, dl, axis = t
                key = f"src{source}_win{window}_B{B}_dl{dl}"
                entry = results.setdefault(key, {})
                entry["axis_" + ("x" if axis == 0 else "y")] = fut.result()
        for key, entry in results.items():
            entry["sym_residual_abs"] = abs(
                entry["axis_x"]["Gamma_BZ_Ang2"]
                + entry["axis_y"]["Gamma_BZ_Ang2"])
            print(key, entry["axis_x"]["Gamma_BZ_Ang2"],
                  entry["axis_y"]["Gamma_BZ_Ang2"],
                  "gap", entry["axis_x"]["min_gap_eV"], flush=True)
    results["_compute_s"] = round(time.time() - t3, 1)
    results["_meta"] = {
        "nk": nk0, "nb": nb0, "kmesh": list(kmesh), "occ": occ,
        "dk_cart_Ang-1": list(map(float, dk_cart)),
        "ABZ_Ang-2": float(ABZ), "A2d_Ang2": float(A2d),
        "weight_Ang-2": float(weight), "muB_eV_T": MU_B,
        "spn_hermitian_dev": float(herm), "spn_max_Sz_diag": float(szd),
        "read_eig_s": round(t1 - t0, 2), "read_spn_s": round(t2 - t1, 2),
        "unit_note": ("Gamma_BZ_Ang2 in 1/(Ang^2 eV T). SI conversion: "
                       "x * (1e-9 m/Ang)^2 / (1.602176634e-19 J/eV) "
                       "= m^2/(J T). Sym residual is |ax_x + ax_y|."),
    }
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--windows", default="36,48,64,80,96")
    ap.add_argument("--Blist", default="0.125,0.25")
    ap.add_argument("--dllist", default="1e-4,2e-4,4e-4")
    ap.add_argument("--source", default="x,y")
    ap.add_argument("--occ", type=int, default=36)
    ap.add_argument("--kmesh", default="21,21")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", default="stage8m_response.json")
    args = ap.parse_args()
    res = run_response(
        args.dir,
        windows=[int(w) for w in args.windows.split(",")],
        blist=[float(b) for b in args.Blist.split(",")],
        dllist=[float(d) for d in args.dllist.split(",")],
        sources=[s for s in args.source.split(",")],
        occ=args.occ,
        kmesh=tuple(int(k) for k in args.kmesh.split(",")),
        jobs=args.jobs)
    with open(args.out, "w") as f:
        json.dump(res, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
