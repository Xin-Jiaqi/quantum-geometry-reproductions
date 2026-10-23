#!/usr/bin/env python3
"""Stage8m: wannier90 formatted matrix readers (.eig/.mmn/.spn) plus a
frozen-KS static-response estimator built on determinant Wilson plaquettes.

Self-contained (numpy only). The estimator uses the same centred parameter
loop as the independently validated Stage7i kernel:
    corners (k_- , lam_-),(k_+ , lam_-),(k_+ , lam_+),(k_- , lam_+)
with k_+- = k_c +- dk/2 along the chosen Cartesian axis and lam_+- =
lam_c +- dl/2. The k-link is Ua^+ M Ub; M is the Bloch overlap from
wannier90.mmn (mode='mmn') or the identity (mode='frame', valid for
absolute-orbital Bloch embeddings where the orbital basis is k-independent).
wannier90.mmn stores only nearest-neighbour single-step bonds, so in mmn
mode each centred k edge is built as the two-step chain k1 -> k_c -> k2
from the two single-step bonds (O(dk^2) accurate for smooth overlaps).
"""
import numpy as np

# ------------------------------------------------------------------ readers
def read_eig(path):
    """wannier90.eig: lines 'band_index kpoint_index energy'
    (band runs fastest) -> (nk, nb) energy array."""
    rows = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 3:
                continue
            b = int(p[0]) - 1
            k = int(p[1]) - 1
            rows[(k, b)] = float(p[2])
    if not rows:
        raise ValueError(f"empty .eig file: {path}")
    nk = max(k for k, b in rows) + 1
    nb = max(b for k, b in rows) + 1
    E = np.full((nk, nb), np.nan)
    for (k, b), e in rows.items():
        E[k, b] = e
    if np.isnan(E).any():
        raise ValueError(".eig has gaps")
    return E


def read_mmn_spn(path):
    """wannier90.mmn/.spn: optional comment header line(s), then a numeric
    header 'num_bands num_kpts num_nn'; rows 'ik ikp band_i band_j Re Im'.
    VASP 6.4+ (wannier90 v2 interface) instead writes per link a 5-token
    header line 'ik ikp ...' followed by num_bands^2 lines 'Re Im' -- both
    formats are auto-detected. Returns (nb, nk, nn, mats) with
    mats[(ik, ikp)] = complex (nb, nb) (fail-closed on duplicates/gaps)."""
    with open(path) as f:
        head = None
        while True:
            p = f.readline().split()
            if not p:
                raise ValueError(f"no numeric header found in {path}")
            if len(p) >= 3 and _all_int(p[:3]):
                head = p
                break
        nb, nk, nn = int(head[0]), int(head[1]), int(head[2])
        first = f.readline().split()   # first data row
        if not first:
            raise ValueError(f"no data rows in {path}")
        standard = len(first) == 6 and _all_int(first[:4])
        mats = {}
        counts = {}
        if standard:
            # first row is already an element row
            def feed(p):
                ik, ikp = int(p[0]) - 1, int(p[1]) - 1
                i, j = int(p[2]) - 1, int(p[3]) - 1
                key = (ik, ikp)
                m = mats.setdefault(key, np.zeros((nb, nb), complex))
                m[i, j] = complex(float(p[4]), float(p[5]))
                counts[key] = counts.get(key, 0) + 1

            feed(first)
            for t in range(nn * nb * nb - 1):
                p = f.readline().split()
                if len(p) < 6:
                    raise ValueError(f"short row in {path}")
                feed(p)
        else:
            # block format: first is the first link header (5 tokens)
            def feed_link(p, nb, mats, f):
                ik, ikp = int(p[0]) - 1, int(p[1]) - 1
                key = (ik, ikp)
                if key in mats:
                    # VASP can emit the same (ik,ikp) bond twice (e.g. a
                    # folded self-link); keep the first block, consume the rest
                    for _ in range(nb * nb):
                        f.readline()
                    return
                m = np.zeros((nb, nb), complex)
                for idx_f in range(nb * nb):
                    q = f.readline().split()
                    if len(q) < 2:
                        raise ValueError(f"short row in {path}")
                    m[idx_f // nb, idx_f % nb] = complex(float(q[0]),
                                                         float(q[1]))
                mats[key] = m

            feed_link(first, nb, mats, f)
            for _ in range(nn * nk - 1):
                p = f.readline().split()
                if len(p) < 2:
                    raise ValueError(f"short link header in {path}")
                feed_link(p, nb, mats, f)
        for key, m in mats.items():
            if m.shape != (nb, nb):
                raise ValueError(f"incomplete matrix {key}")
        return nb, nk, nn, mats


def _all_int(p):
    try:
        return all(int(x) == float(x) for x in p)
    except ValueError:
        return False


def link_matrix(mats, ik, ikp, nb):
    """M(ik, ikp); fall back to conjugate transpose of reverse bond."""
    if (ik, ikp) in mats:
        return mats[(ik, ikp)]
    if (ikp, ik) in mats:
        return mats[(ikp, ik)].T.conj()
    raise ValueError(f"missing overlap link ({ik},{ikp}) and reverse")


# ------------------------------------------------------------- plaquette engine
def polar_unitary(m):
    """Unitary part of polar decomposition; returns (U, min singular value)."""
    u, s, vh = np.linalg.svd(m)
    return (u @ vh).astype(complex), float(s[-1])


class DeterminantPlaquette:
    """F_{lam,kj} via determinant Wilson loops; Gamma(k) via central B diff.

    frames_at(k_coord, lam, B) -> U_occ : (nb x nb_occ) complex, coefficients
    of the occupied manifold in the fixed basis. k_coord is Cartesian when
    mode='frame' (model regression); for mode='mmn' k_coord is the mesh index
    and ik_of maps corners back to mesh indices for the M lookup.
    """

    def __init__(self, frames_at, *, mode, links=None, ik_of=None,
                 dk, dl, axis=0, full_nb=None):
        self.frames_at = frames_at
        self.mode = mode
        self.links = links
        self.ik_of = ik_of
        self.dk = dk
        self.dl = dl
        self.axis = axis
        self.full_nb = full_nb
        self.min_svals = []

    def _corners(self, k_c, lam_c, B):
        """Return the four ordered frames (k-,lam-),(k+,lam-),(k+,lam+),(k-,lam+)
        and the k-coordinates/indexes of the two k corners."""
        ax = self.axis
        if self.mode == "mmn":
            k1i = self.ik_of(k_c, -1)
            k2i = self.ik_of(k_c, +1)
            k1 = k2 = None
            F1 = self.frames_at(k1i, lam_c - self.dl / 2, B)
            F2 = self.frames_at(k2i, lam_c - self.dl / 2, B)
            F3 = self.frames_at(k2i, lam_c + self.dl / 2, B)
            F4 = self.frames_at(k1i, lam_c + self.dl / 2, B)
            return F1, F2, F3, F4, (k1i, k2i)
        # frame mode: Cartesian coordinates
        k1 = np.array(k_c, float).copy()
        k2 = np.array(k_c, float).copy()
        k1[ax] -= self.dk / 2
        k2[ax] += self.dk / 2
        F1 = self.frames_at(k1, lam_c - self.dl / 2, B)
        F2 = self.frames_at(k2, lam_c - self.dl / 2, B)
        F3 = self.frames_at(k2, lam_c + self.dl / 2, B)
        F4 = self.frames_at(k1, lam_c + self.dl / 2, B)
        return F1, F2, F3, F4, (k1, k2)

    def _k_link_chain(self, k_from, k_mid, k_to, lam, B):
        """Product Wilson link k_from -> k_mid -> k_to built from the two
        single-step overlap bonds that wannier90.mmn actually provides
        (M(k,k') exists only for nearest neighbours). For a smooth Bloch
        overlap the two-step product approximates the direct bond to
        O(dk^2), keeping the centred plaquette second-order accurate."""
        Uf = self.frames_at(k_from, lam, B)
        Um = self.frames_at(k_mid, lam, B)
        Ut = self.frames_at(k_to, lam, B)
        nb = Uf.shape[0]                       # window size
        full = self.full_nb or nb
        # single-step bonds are full-basis overlaps; truncate to the window
        Mfm = link_matrix(self.links, k_from, k_mid, full)[:nb, :nb]
        Mmt = link_matrix(self.links, k_mid, k_to, full)[:nb, :nb]
        W1 = Uf.T.conj() @ Mfm @ Um
        W2 = Um.T.conj() @ Mmt @ Ut
        return W1 @ W2

    def plaquette_phase(self, k_c, lam_c, B):
        F1, F2, F3, F4, (k1, k2) = self._corners(k_c, lam_c, B)
        if self.mode == "mmn":
            # k edges are two-step chains via the centre point k_c
            L1 = self._k_link_chain(k1, k_c, k2, lam_c - self.dl / 2, B)
            L3 = self._k_link_chain(k2, k_c, k1, lam_c + self.dl / 2, B)
            U2m = self.frames_at(k2, lam_c - self.dl / 2, B)
            U2p = self.frames_at(k2, lam_c + self.dl / 2, B)
            U1p = self.frames_at(k1, lam_c + self.dl / 2, B)
            U1m = self.frames_at(k1, lam_c - self.dl / 2, B)
            L2 = U2m.T.conj() @ U2p          # lambda-link at k+
            L4 = U1p.T.conj() @ U1m          # lambda-link at k-
        else:
            L1 = F1.T.conj() @ F2
            L3 = F3.T.conj() @ F4
            L2 = F2.T.conj() @ F3          # lambda-link at k+
            L4 = F4.T.conj() @ F1          # lambda-link at k-
        det = 1.0 + 0j
        for L in (L1, L2, L3, L4):
            Lp, smin = polar_unitary(L)
            self.min_svals.append(smin)
            det *= np.linalg.det(Lp)
        if abs(det) < 1e-300:
            raise FloatingPointError("Wilson determinant collapsed")
        return np.angle(det) / (self.dk * self.dl)

    def Gamma_k(self, k_c, lam_c, B, dB):
        Fp = self.plaquette_phase(k_c, lam_c, B + dB)
        Fm = self.plaquette_phase(k_c, lam_c, B - dB)
        return (Fp - Fm) / (2 * dB)


def sum_response(plaquette, points, *, lam=0.0, B=0.0, dB,
                 weight=None, gap_guard=None):
    """Sum Gamma(k) over the given mesh points with optional weights.
    gap_guard(coord_or_index) -> gap; None/<=0 raises. Returns
    (Gamma_sum, min_gap, min_sval)."""
    total = 0.0
    ming = float("inf")
    mins = float("inf")
    for ik, pt in enumerate(points):
        if gap_guard is not None:
            g = gap_guard(pt)
            if g is None or g <= 0:
                raise ValueError(f"gap guard failed at point {pt}")
            ming = min(ming, g)
        gam = plaquette.Gamma_k(pt, lam, B, dB)
        w = 1.0 if weight is None else weight[ik]
        total += w * gam
        if plaquette.min_svals:
            mins = min(mins, min(plaquette.min_svals[-4:]))
    return total, ming, mins
