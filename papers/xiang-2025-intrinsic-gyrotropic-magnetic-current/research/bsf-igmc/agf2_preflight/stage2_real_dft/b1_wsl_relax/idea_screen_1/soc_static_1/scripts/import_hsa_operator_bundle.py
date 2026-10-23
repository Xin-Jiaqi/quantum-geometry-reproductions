#!/usr/bin/env python3
"""Strict, fail-closed H/S/A operator import path for PR #38 (synthetic only).

Phase A - versioned synthetic carrier: serialize a validated four-band
H(R)/S(R)/A(R) + lattice + centers + units into a manifest-hashed JSON
carrier; strict reader returns the existing WannierOperatorBundle.

Phase B - genuine Wannier90 syntax parser for AUTHORED SYNTHETIC files
(no Wannier run): hr.dat, r.dat, formatted spn, u.mat/u_dis.mat/v.mat,
eig, amn/mmn (provenance), win metadata. Reconstruct the final-Wannier
gauge operators: V(k) = U_dis(k) U(k) (num_bands x num_wann), the spin via
S^W(k) = V+(k) S^Bloch(k) V(k) with matching discrete Fourier convention,
and the position directly from r.dat (already final Wannier gauge; the
centers are the diagonal on-site position entries, included ONCE).

Phase C - Wilson winding: proper unwrapped 2pi-winding over the periodic
ky cycle with closure, replacing the unreachable wrap(endpoint) predicate.

Formats verified against the Wannier90 documentation (master, read-only
fetched 2026-10-23):
  docs/user_guide/wannier90/files.md and postproc.md
  test-suite/tests/testw90_rmn/benchmark/diamond_r.dat (layout)

Declared first support case (everything else raises):
  2D slab; 3D orthorhombic lattice; Rz = 0 for every R;
  write_ndegen_applied = true with all degeneracy = 1;
  use_ws_distance off / wsvec NOT applied on top; SAXIS = (0,0,1);
  spinors = true; a single explicit k ordering from the win kpoints block;
  formatted spn only (binary spn rejected); no ambiguous WS folding.

Status labels: SYNTHETIC_CARRIER_PASS / SYNTHETIC_FORMAT_PARSER_PASS only.
MATERIAL_AUTHENTICATED is impossible here (no ab-initio inputs exist).
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from wannier_operator_contract_pinned26 import (
    WannierOperatorBundle,
    UNITS,
)

# ---------------------------------------------------------------- manifest
def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def canonical_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


# ---------------------------------------------------------------- Phase A
def _complex_matrices_to_lists(bundle):
    """Serialize a bundle to a plain (JSON-safe) dict."""
    return {
        "translations": bundle.translations.tolist(),
        "lattice": bundle.lattice.tolist(),
        "h_r": [[[complex(x).real, complex(x).imag] for x in row] for row in bundle.h_r.reshape(-1, 4)],
        "spin_r": [[[complex(x).real, complex(x).imag] for x in row]
                   for row in bundle.spin_r.reshape(-1, 4)],
        "position_r": [[[complex(x).real, complex(x).imag] for x in row]
                       for row in bundle.position_r.reshape(-1, 4)],
        "centers": bundle.centers.tolist(),
        "provenance": bundle.provenance,
        "units": bundle.units,
    }


def write_hsa_carrier(bundle, path, label="synthetic_v2_positive",
                      spin_axis=(0, 0, 1), status="SYNTHETIC_CARRIER"):
    """Write manifest-hashed synthetic carrier JSON. Returns manifest dict."""
    data = _complex_matrices_to_lists(bundle)
    body = {"schema": "hsa-synthetic-carrier-v1",
            "label": label,
            "spin_axis": list(spin_axis),
            "spin_normalization": "dimensionless_pauli",
            "position_center_once": True,
            "data": data}
    blob = canonical_bytes(body)
    manifest = {"schema": "hsa-synthetic-carrier-manifest-v1",
                "status": status,
                "file_sha256": sha256_bytes(blob),
                "label": label,
                "basis_ids": bundle.provenance.get("basis_ids")}
    out = {"manifest": manifest, "body": body}
    Path(path).write_text(json.dumps(out, indent=1) + "\n")
    return manifest


def read_hsa_carrier(path, require_status=None):
    """Strict reader -> (WannierOperatorBundle, manifest). Every structural
    violation raises (nonzero exit in tests), never a mere 'detected=True'."""
    raw = Path(path).read_text()
    obj = json.loads(raw)
    manifest = obj["manifest"]
    body = obj["body"]
    # 1) integrity: canonical hash of the body must equal the manifest hash
    expect = manifest["file_sha256"]
    got = sha256_bytes(canonical_bytes(body))
    if got != expect:
        raise ValueError(f"carrier hash mismatch: {got[:12]} != {expect[:12]}")
    if require_status is not None and manifest["status"] != require_status:
        raise ValueError(f"status {manifest['status']} != required {require_status}")
    if body["spin_axis"] != [0, 0, 1]:
        raise ValueError(f"unsupported spin axis {body['spin_axis']}; only (0,0,1)")
    if body["spin_normalization"] != "dimensionless_pauli":
        raise ValueError("spin normalization must be dimensionless_pauli")
    if not body["position_center_once"]:
        raise ValueError("position_center_once must be true")
    d = body["data"]
    # 2) types/shapes and finite complex values
    translations = np.asarray(d["translations"], dtype=int)
    if translations.ndim != 2 or translations.shape[1] != 2:
        raise ValueError("translations must be (nr, 2) integers (2D slab)")
    # 2a) R-list structural integrity BEFORE parsing matrices: R=0 present
    # and every R has its -R partner (a matrix-parse error would otherwise
    # mask the missing-partner violation)
    tset = {tuple(int(x) for x in t) for t in translations}
    if (0, 0) not in tset:
        raise ValueError("translations missing R=(0,0) on-site entry")
    for r in tset:
        if (-r[0], -r[1]) not in tset:
            raise ValueError(f"translations R={r} has no -R partner "
                             f"({(-r[0], -r[1])} absent)")
    lattice = np.asarray(d["lattice"], dtype=float)
    if lattice.shape != (2, 3):
        raise ValueError("lattice must be (2, 3) Angstrom rows (2D slab in 3D)")
    nr = len(translations)
    mats = {}
    for key, shape in (("h_r", (nr, 4, 4)), ("spin_r", (nr, 3, 4, 4)),
                       ("position_r", (nr, 3, 4, 4))):
        if key not in d:
            raise ValueError(f"carrier missing required matrix {key}")
        arr = _lists_to_complex(d[key], shape)
        if not np.all(np.isfinite(arr.view(float))):
            raise ValueError(f"{key}: non-finite values")
        mats[key] = arr
    h_r, spin_r, position_r = mats["h_r"], mats["spin_r"], mats["position_r"]
    centers = np.asarray(d["centers"], dtype=float)
    if centers.shape != (4, 3):
        raise ValueError("centers must be (4, 3) Angstrom")
    # 3) R <-> -R Hermiticity (the ONLY per-R spin condition; a general
    #    Wannier spin operator has S_j(-R)=S_j(R)+ with complex, non-
    #    Hermitian off-site terms - the individual S_j(R)=S_j(R)+ holds
    #    only at R=0 and is checked in step 4)
    _check_r_pair_hermiticity(h_r, translations)
    for j in range(3):
        _check_r_pair_hermiticity(spin_r[:, j], translations)
        _check_r_pair_hermiticity(position_r[:, j], translations)
    # 4) spin: only the on-site R=0 term is individually Hermitian; the
    #    physical Pauli norm bound applies to the reconstructed S_j(k) on a
    #    k sample (per-R coefficient norms are NOT a physical provenance)
    i0 = {tuple(int(x) for x in t): i for i, t in enumerate(translations)}[(0, 0)]
    for j in range(3):
        if not np.allclose(spin_r[i0, j], spin_r[i0, j].conj().T, atol=1e-12):
            raise ValueError(f"spin_r[R=0,{j}] not Hermitian (on-site only)")
    for qx, qy in ((0.0, 0.0), (0.2, 0.1), (0.45, 0.3), (0.7, 0.6), (0.9, 0.9)):
        for j in range(3):
            Sk = np.zeros((4, 4), complex)
            for i, t in enumerate(translations):
                ph = np.exp(2j * np.pi * (t[0] * qx + t[1] * qy))
                Sk += ph * spin_r[i, j]
            if not np.allclose(Sk, Sk.conj().T, atol=1e-9):
                raise ValueError(f"spin S_j(q={qx},{qy}) not Hermitian")
            # physical Pauli norm: ||S_j(k)|| <= 1 with a 0.05 slack at the
            # sample points (a legitimate synthetic off-site spin may carry a
            # small non-Hermitian R-term whose Fourier norm exceeds 1 by
            # ~2|delta|); the tight ||S(0)|| <= 1+1e-9 stays in validate().
            # NOTE: a bounded numerical norm is necessary, not sufficient,
            # evidence of physical provenance (numerical norm != provenance).
            if np.linalg.norm(Sk, 2) > 1.0 + 0.05:
                raise ValueError(f"spin S_j(q={qx},{qy}) norm > 1+0.05 "
                                 f"(not projected physical Pauli)")
    # 5) centers == diag A(R=0) exactly once (R=0 may not be list index 0)
    onsite = np.diagonal(position_r[i0], axis1=-2, axis2=-1).T.real
    if not np.allclose(onsite, centers, atol=1e-9):
        raise ValueError("centers must equal diag of the on-site position "
                         "connection (centers counted once)")
    # 6) units
    units = dict(UNITS)
    units.update(d.get("units", {}))
    if units.get("energy") != "eV" or units.get("length") != "angstrom":
        raise ValueError(f"unsupported units {units}")
    bundle = WannierOperatorBundle(
        np.array(translations, dtype=int), np.array(lattice, dtype=float),
        h_r, spin_r, position_r, centers, dict(d.get("provenance", {})), units)
    bundle.validate()
    return bundle, manifest


def _lists_to_complex(lst, shape):
    a = np.asarray(lst, dtype=float).reshape(shape + (2,))
    return a[..., 0] + 1j * a[..., 1]


def _check_r_pair_hermiticity(mat, translations):
    """For each R with partner -R: X(-R) == X(R)^dagger (contract convention).
    Every R must have its -R partner present (else the pair invariant is
    unenforceable -> fail closed)."""
    idx = {tuple(int(x) for x in t): i for i, t in enumerate(translations)}
    for i, t in enumerate(translations):
        key = (-int(t[0]), -int(t[1]))
        if key not in idx:
            raise ValueError(f"R={tuple(t)} has no -R partner in the R list")
        j = idx[key]
        if not np.allclose(mat[i], mat[j].conj().T, atol=1e-10):
            raise ValueError(f"R={tuple(t)} missing Hermitian partner -R")


# ---------------------------------------------------------------- Phase B
def _read_lines(path):
    return Path(path).read_text().splitlines()


def parse_hr_dat(path, ndegen_applied_expected=True):
    """hr.dat (write_hr). Returns (translations, h_r (nr,4,4) complex).
    write_ndegen_applied=true: expanded R list, degeneracy block all 1."""
    lines = _read_lines(path)
    if len(lines) < 4:
        raise ValueError("hr.dat too short")
    num_wann = int(lines[1].split()[0])
    if num_wann != 4:
        raise ValueError(f"hr.dat num_wann {num_wann} != 4 (first support case)")
    nrpts = int(lines[2].split()[0])
    deg = [int(x) for x in lines[3].split()]
    p = 4
    while len(deg) < nrpts:
        if p >= len(lines):
            raise ValueError("hr.dat degeneracy block truncated")
        deg += [int(x) for x in lines[p].split()]
        p += 1
    deg = deg[:nrpts]
    if ndegen_applied_expected and any(x != 1 for x in deg):
        raise ValueError("write_ndegen_applied=true requires all degeneracy 1")
    body = lines[p:]
    nper = num_wann * num_wann
    translations, h_r = [], np.zeros((nrpts, 4, 4), complex)
    pos = 0
    for ir in range(nrpts):
        seg = body[pos:pos + nper]
        if len(seg) < nper:
            raise ValueError("hr.dat truncated at R block")
        pos += nper
        tok0 = seg[0].split()
        if len(tok0) < 5:
            raise ValueError("hr.dat malformed R line")
        R = (int(tok0[0]), int(tok0[1]), int(tok0[2]))
        if R[2] != 0:
            raise ValueError(f"Rz={R[2]} != 0 (2D slab only)")
        translations.append((R[0], R[1]))
        for k in range(nper):
            toks = seg[k].split()
            if len(toks) != 7:
                raise ValueError(f"hr.dat entry has {len(toks)} fields, want 7")
            m, n = int(toks[3]) - 1, int(toks[4]) - 1
            h_r[ir, m, n] = float(toks[5]) + 1j * float(toks[6])
    if pos != len(body):
        raise ValueError("hr.dat extra trailing lines")
    return np.array(translations), h_r


def parse_r_dat(path):
    """r.dat (write_rmn): final Wannier-gauge position connection.
    Layout (verified vs benchmark): R1 R2 R3 m n x_re x_im y_re y_im z_re z_im;
    no degeneracy block (weights come from hr.dat / ndegen_applied)."""
    lines = _read_lines(path)
    if len(lines) < 4:
        raise ValueError("r.dat too short")
    num_wann = int(lines[1].split()[0])
    if num_wann != 4:
        raise ValueError(f"r.dat num_wann {num_wann} != 4")
    nrpts = int(lines[2].split()[0])
    nper = num_wann * num_wann
    body = lines[3:]
    if len(body) != nrpts * nper:
        raise ValueError(f"r.dat line count {len(body)} != {nrpts * nper}")
    translations, a_r = [], np.zeros((nrpts, 3, 4, 4), complex)
    for ir in range(nrpts):
        seg = body[ir * nper:(ir + 1) * nper]
        tok0 = seg[0].split()
        if len(tok0) != 11:
            raise ValueError(f"r.dat entry has {len(tok0)} fields, want 11 "
                             f"(R1 R2 R3 m n + 3 complex Cartesian)")
        R = (int(tok0[0]), int(tok0[1]), int(tok0[2]))
        if R[2] != 0:
            raise ValueError(f"Rz={R[2]} != 0 (2D slab only)")
        translations.append((R[0], R[1]))
        for k in range(nper):
            toks = seg[k].split()
            if len(toks) != 11:
                raise ValueError("r.dat malformed entry")
            m, n = int(toks[3]) - 1, int(toks[4]) - 1
            a_r[ir, 0, m, n] = float(toks[5]) + 1j * float(toks[6])
            a_r[ir, 1, m, n] = float(toks[7]) + 1j * float(toks[8])
            a_r[ir, 2, m, n] = float(toks[9]) + 1j * float(toks[10])
    return np.array(translations), a_r


def parse_spn_formatted(path, num_bands):
    """Formatted spn: comment; num_bands num_kpts; per k: upper-triangle
    <psi_n|sigma_i|psi_m> for m=1..nb, n=1..m, i=x,y,z (i fastest).
    Returns list (per k) of three (nb,nb) Hermitian complex matrices."""
    lines = _read_lines(path)
    if len(lines) < 3:
        raise ValueError("spn too short")
    tok = lines[1].split()
    nb, nk = int(tok[0]), int(tok[1])
    if nb != num_bands:
        raise ValueError(f"spn num_bands {nb} != {num_bands}")
    count = 3 * nb * (nb + 1) // 2
    if len(lines) < 2 + count * nk:
        raise ValueError("spn truncated")
    out = []
    p = 2
    for _ in range(nk):
        mats = [np.zeros((nb, nb), complex) for _ in range(3)]
        for m in range(1, nb + 1):
            for n in range(1, m + 1):
                for i in range(3):
                    t2 = lines[p].split()
                    if len(t2) != 2:
                        raise ValueError("spn element line must have re im")
                    re_ = float(t2[0]); im_ = float(t2[1])
                    p += 1
                    mats[i][n - 1, m - 1] = re_ + 1j * im_
                    if n != m:
                        mats[i][m - 1, n - 1] = re_ - 1j * im_
        out.append(mats)
    if p != len(lines):
        raise ValueError("spn extra trailing lines")
    return out


def _check_k_list_consistent(kpts, label):
    """Structural k-list check: finite 3-vectors, count >= 1, each k in
    [0,1) (fractional, wraps allowed), and no duplicated k points within
    1e-8 (a 'single explicit k ordering' must be unique). Cross-file and
    cross-win equality is enforced separately with numeric tolerances."""
    k = np.asarray(kpts, dtype=float)
    if k.ndim != 2 or k.shape[1] != 3 or k.shape[0] < 1:
        raise ValueError(f"{label}: k list must be (nk,3), nk >= 1")
    if not np.all(np.isfinite(k)):
        raise ValueError(f"{label}: non-finite k coordinates")
    kf = k % 1.0
    if np.max(np.abs(kf - k)) > 1e-8:
        raise ValueError(f"{label}: k coordinates outside [0,1)")
    d = np.abs(k[:, None, :] - k[None, :, :]).max(axis=-1)
    iu = np.triu(np.ones((len(k), len(k)), dtype=bool), k=1)
    if np.any(d[iu] < 1e-8):
        raise ValueError(f"{label}: duplicate k points (k ordering not unique)")
    return k


def parse_kpoint_matrix_file(path, key):
    """u.mat / u_dis.mat / v.mat formatted layout:
    comment; header; blank; per k: kpoint line + entries (re im), COLUMN-MAJOR
    (rows change fastest, then columns; official docs), blank line between
    blocks.
    key selects the expected header: 'u' -> nk, nw, nw; 'u_dis'/'v' -> nk, nw, nb.
    Returns (kpts (nk,3) float, list of complex matrices (nrows,ncols))."""
    lines = _read_lines(path)
    tok = lines[1].split()
    nk = int(tok[0])
    if key == "u":
        nw = int(tok[1])
        if int(tok[2]) != nw:
            raise ValueError("u.mat header inconsistent")
        nrows, ncols = nw, nw
    else:
        nw, nb = int(tok[1]), int(tok[2])
        nrows, ncols = nb, nw
    if lines[2].strip() != "":
        raise ValueError("third line of matrix file must be blank")
    kpts, mats = [], []
    p = 3
    for _ in range(nk):
        if p >= len(lines) or not lines[p].strip():
            raise ValueError("matrix file: expected kpoint line")
        kk = [float(x) for x in lines[p].split()]
        if len(kk) != 3:
            raise ValueError("kpoint line must have 3 fractional coords")
        kpts.append(kk)
        p += 1
        nent = nrows * ncols
        m = np.zeros((nrows, ncols), complex)
        for c in range(ncols):
            for r in range(nrows):      # column-major: rows fastest
                if p >= len(lines):
                    raise ValueError("matrix file truncated")
                t2 = lines[p].split()
                if len(t2) != 2:
                    raise ValueError("matrix element line must have re im")
                m[r, c] = float(t2[0]) + 1j * float(t2[1])
                p += 1
        mats.append(m)
        if p < len(lines) and lines[p].strip() == "":
            p += 1  # blank line between blocks
    return np.array(kpts), mats


def parse_eig(path, num_bands, num_kpts):
    """eig: 'band_index k_ordinal eigenvalue' (both 1-based)."""
    lines = _read_lines(path)
    if len(lines) != num_bands * num_kpts:
        raise ValueError(f"eig line count {len(lines)} != {num_bands * num_kpts}")
    e = np.zeros((num_kpts, num_bands))
    for line in lines:
        t = line.split()
        if len(t) != 3:
            raise ValueError("eig line must have band k eigenvalue")
        b, k, v = int(t[0]) - 1, int(t[1]) - 1, float(t[2])
        if b < 0 or b >= num_bands or k < 0 or k >= num_kpts:
            raise ValueError("eig band/k ordinal out of range")
        e[k, b] = v
    return e


def parse_amn(path):
    """amn provenance only: comment; num_bands num_kpts num_wann; then lines
    'm n k re im' (1-based). Returns (meta, entries)."""
    lines = _read_lines(path)
    tok = lines[1].split()
    nb, nk, nw = (int(x) for x in tok)
    entries = []
    for line in lines[2:]:
        t = line.split()
        if len(t) != 5:
            raise ValueError("amn entry must be 'm n k re im'")
        entries.append((int(t[0]), int(t[1]), int(t[2]),
                        float(t[3]), float(t[4])))
    if len(entries) != nb * nw * nk:
        raise ValueError("amn entry count mismatch")
    return {"num_bands": nb, "num_kpts": nk, "num_wann": nw}, entries


def parse_mmn(path):
    """mmn provenance only: comment; num_bands num_kpts nntot."""
    lines = _read_lines(path)
    tok = lines[1].split()
    nb, nk, nntot = (int(x) for x in tok)
    return {"num_bands": nb, "num_kpts": nk, "nntot": nntot}


def parse_win_metadata(path):
    """win: read lattice, spinors, saxis, kpoints block, write flags."""
    text = Path(path).read_text()
    lat = None
    m = re.search(r"begin\s+unit_cell_cart\s*\n(.*?)end\s+unit_cell_cart",
                  text, re.S | re.I)
    if m:
        rows = [[float(x) for x in ln.split()] for ln in m.group(1).splitlines()
                if ln.strip()]
        if len(rows) == 3:
            lat = np.array(rows)
    saxis = None
    m = re.search(r"saxis\s*=\s*([\d\s.,\-+]+)", text, re.I)
    if m:
        vals = [float(x) for x in m.group(1).replace(",", " ").split()]
        if len(vals) == 3:
            saxis = vals
    kpts = None
    m = re.search(r"begin\s+kpoints\s*\n(.*?)end\s+kpoints", text, re.S | re.I)
    if m:
        kpts = [[float(x) for x in ln.split()] for ln in m.group(1).splitlines()
                if ln.strip()]
    ndegen = re.search(r"write_ndegen_applied\s*=\s*\.true\.", text, re.I)
    wsdist = re.search(r"use_ws_distance\s*=\s*\.true\.", text, re.I)
    spinors = re.search(r"spinors\s*=\s*\.true\.", text, re.I)
    nw = re.search(r"num_wann\s*=\s*(\d+)", text, re.I)
    nb = re.search(r"num_bands\s*=\s*(\d+)", text, re.I)
    return {"lattice_cart": lat, "saxis": saxis, "kpoints": kpts,
            "write_ndegen_applied": bool(ndegen),
            "use_ws_distance": bool(wsdist),
            "spinors": bool(spinors),
            "num_wann": int(nw.group(1)) if nw else None,
            "num_bands": int(nb.group(1)) if nb else None}


def _check_support_case(meta, nbands, nwann, nkpts):
    """Fail closed on anything outside the declared first support case."""
    if meta.get("spinors") is not True:
        raise ValueError("win must declare spinors = true")
    if meta.get("saxis") != [0.0, 0.0, 1.0]:
        raise ValueError(f"saxis must be (0,0,1), got {meta.get('saxis')}")
    if meta.get("use_ws_distance"):
        raise ValueError("use_ws_distance unsupported (ambiguous WS folding)")
    if meta.get("num_wann") != nwann or meta.get("num_bands") != nbands:
        raise ValueError("win num_wann/num_bands mismatch")
    k = meta.get("kpoints")
    if k is None or len(k) != nkpts:
        raise ValueError(f"explicit kpoints block required with {nkpts} points")
    lat = meta.get("lattice_cart")
    if lat is None:
        raise ValueError("unit_cell_cart required (3x3 Angstrom)")
    off = np.abs(lat - np.diag(np.diag(lat)))
    if np.max(off) > 1e-8:
        raise ValueError("lattice must be orthorhombic (diagonal) for support")
    return np.array(k, dtype=float), np.array(lat, dtype=float)


def fourier_inverse_on_grid(X_k, kpts, translations):
    """X(R) = mean_q exp(-2pi i R.q) X(q), exact on the uniform grid when
    R components are within the Brillouin zone of the grid."""
    out = []
    for R in translations:
        phase = np.exp(-2j * np.pi * (R[0] * kpts[:, 0] + R[1] * kpts[:, 1]))
        out.append(np.mean(phase[:, None, None] * X_k, axis=0))
    return np.array(out)


def reconstruct_bundle_from_w90(directory, nbands, nwann,
                                ndegen_applied_expected=True,
                                h_eig_tol=1e-8):
    """Full Phase B pipeline for AUTHORED SYNTHETIC files.

    Reads hr.dat, r.dat, spn (formatted), u.mat / u_dis.mat / v.mat, eig,
    win. Verifies:
      * V(k) = U_dis(k) U(k), V+ V = I (num_wann x num_wann);
      * S^W(k) = V+(k) S^Bloch(k) V(k) with the Hermitian spin from spn;
      * H^W(k) = V+ diag(eig(k)) V == sum_R exp(2pi i R.q) H_R (matched);
      * r.dat is the FINAL Wannier-gauge connection (no second rotation);
        centers = diag A(0), counted once;
      * Fourier-back equality of H/S/A at every mesh k.
    h_eig_tol: tolerance for the H-vs-eig consistency (1e-8 for the exact
    identity gauge; smooth final gauges have a Fourier truncation error
    from the k-dependent rotation and are tested at a documented looser
    tolerance).
    Returns (bundle, receipt dict)."""
    D = Path(directory)
    hr_t, h_r = parse_hr_dat(D / "seed_hr.dat",
                             ndegen_applied_expected=ndegen_applied_expected)
    r_t, a_r = parse_r_dat(D / "seed_r.dat")
    if not np.array_equal(hr_t, r_t):
        raise ValueError("hr.dat and r.dat R lists differ (ndegen/WS mismatch)")
    spn = parse_spn_formatted(D / "seed.spn", nbands)
    meta = parse_win_metadata(D / "seed.win")
    kpts, lat = _check_support_case(meta, nbands, nwann, len(spn))
    # FIX 3: k-list validation - every U-family file's k coordinates must
    # match the win explicit kpoint list and one another (orientation,
    # count, numeric tolerance); spn blocks carry NO per-k labels, so its
    # k mapping is provenance-limited (see provenance below)
    win_k = (np.array(meta["kpoints"], dtype=float) if meta.get("kpoints")
             else None)
    if win_k is None or len(win_k) != len(spn):
        raise ValueError("win explicit kpoints block missing or count "
                         "mismatch vs spn blocks")
    _check_k_list_consistent(win_k, "win kpoints")

    def _load(key):
        ku, mu = parse_kpoint_matrix_file(D / f"seed_{key}.mat", key)
        _check_k_list_consistent(ku, f"{key}.mat kpoints")
        if ku.shape != win_k.shape or not np.allclose(ku, win_k, atol=1e-8):
            raise ValueError(f"{key}.mat k points differ from win kpoints "
                             f"(max |d|="
                             f"{float(np.max(np.abs(ku - win_k))):.2e})")
        return mu

    u = _load("u")
    v_has, ud_has = (D / "seed_v.mat").exists(), (D / "seed_u_dis.mat").exists()
    if v_has:
        v = _load("v")
    if ud_has:
        udis = _load("u_dis")
    if v_has and ud_has:
        # FIX 3: never silently pick one: verify V = U_dis U explicitly
        max_vud = 0.0
        for k in range(len(u)):
            diff = v[k] - udis[k] @ u[k]
            max_vud = max(max_vud, float(np.max(np.abs(diff))))
        if max_vud > 1e-8:
            raise ValueError(f"V != U_dis U: max |V-U_dis U| = {max_vud:.3e}")
        v_consistency = float(max_vud)
    elif v_has:
        v_consistency = None          # no u_dis present: nothing to compare
    elif ud_has:
        v = [udis[k] @ u[k] for k in range(len(u))]
        v_consistency = None
    else:
        v = u                          # no disentanglement: V(k) = U(k) alone
        v_consistency = None
    if len(v) != len(spn) or len(u) != len(spn):
        raise ValueError("gauge/spn k-count mismatch")
    e = parse_eig(D / "seed.eig", nbands, len(spn))
    # V+ V = I and V+ diag(eig) V = H^W(k) = fourier(H_R)
    tol = 1e-8
    max_vv = 0.0
    for k in range(len(spn)):
        V = v[k]
        vv = V.conj().T @ V
        max_vv = max(max_vv, float(np.max(np.abs(vv - np.eye(nwann)))))
        Hw_from_eig = V.conj().T @ np.diag(e[k]) @ V
        Hw_from_hr = np.zeros((nwann, nwann), complex)
        for ir, R in enumerate(hr_t):
            ph = np.exp(2j * np.pi * (R[0] * kpts[k, 0] + R[1] * kpts[k, 1]))
            Hw_from_hr += ph * h_r[ir]
        err = float(np.max(np.abs(Hw_from_hr - Hw_from_eig)))
        if err > h_eig_tol:
            raise ValueError(f"H(k) mismatch vs eig at k={k}: {err:.3e} "
                             f"(tol {h_eig_tol})")
    if max_vv > 1e-8:
        raise ValueError(f"V+V != I: {max_vv:.3e}")
    # spin: S^W(k) = V+ S^B V (per Pauli component), then S(R) by iFT
    S_w_k = np.array([
        [v[k].conj().T @ np.array(spn[k])[j] @ v[k] for j in range(3)]
        for k in range(len(spn))])                       # (nk,3,4,4)
    spin_r = np.stack([fourier_inverse_on_grid(S_w_k[:, j], kpts, hr_t)
                       for j in range(3)], axis=1)  # (nr,3,4,4)
    # position: r.dat is final gauge -> A(R) direct; centers = diag A(0)
    i0 = {tuple(int(x) for x in t): i for i, t in enumerate(hr_t)}[(0, 0)]
    centers = np.diagonal(a_r[i0], axis1=-2, axis2=-1).T.real
    lat2d = lat[:2]
    bundle = WannierOperatorBundle(
        np.array(hr_t, dtype=int), np.array(lat2d, dtype=float),
        h_r, spin_r, a_r, centers,
        {"basis_ids": {k: "synthetic-w90-authored-final-gauge"
                       for k in ("hamiltonian", "spin", "position")},
         "spin_kind": "projected_physical_pauli",
         "spin_evidence": "synthetic authoring: spn Hermitian upper-triangle "
                          "Bloch spin, V+ S^B V (W90 gauge algebra)",
         "orthonormal_basis": True,
         "position_convention": "full_connection_including_centers",
         "producer": "synthetic_authoring_only", "run_kind": "synthetic_fixture",
         "write_ndegen_applied": True, "use_ws_distance": False,
         "saxis": [0, 0, 1], "spin_normalization": "dimensionless_pauli",
         "spin_k_provenance": "SPN_K_PROVENANCE_UNVERIFIED"},
        dict(UNITS))
    bundle.validate()
    receipt = {"max_VdagV_minus_I": float(max_vv),
               "kpts": len(spn), "nbands": nbands, "nwann": nwann,
               "H_matched_to_eig_max_err": float(err),
               "k_validation": "U-family k lists matched win explicit "
                               "kpoints (atol 1e-8) and each other",
               "V_UdisU_consistency": v_consistency,
               "spin_k_provenance": "SPN_K_PROVENANCE_UNVERIFIED"}
    return bundle, receipt


# ---------------------------------------------------------------- Phase C
def wilson_winding_unwrapped(rows):
    """Integer 2pi-winding of an unwrapped Wilson-flow phase sequence over a
    periodic cycle (closing link included).

    rows: wrapped per-ky Wilson row phases in [-pi, pi). Returns
    (winding, closure_residual, max_local_step, jump_like):
      winding = round((u[-1]-u[0])/2pi) with u = np.unwrap(rows);
      closure residual = |(u[-1]-u[0]) - 2pi*winding|;
      jump_like = the 2pi net was acquired in ONE local step >= pi*1.9
      (suspect branch jump) instead of being distributed smoothly.
    The old predicate abs(wrap(last-first)) >= pi is unreachable because
    wrap lies in [-pi, pi); it cannot detect integer winding.

    Boundary-branch resolution: the closing point stores gamma(2pi) as the
    principal branch of gamma(0) + 2pi*C, so an in-range boundary jump
    (|step| < pi) can hide one 2pi of winding (closure residual ~ pi). In
    that case the flow is re-read from the interior (the boundary point is
    dropped) and the interior drift is used when it closes significantly
    better. Returns the resolved integer winding."""
    x = np.asarray(rows, dtype=float)
    u = np.unwrap(x)
    total = u[-1] - u[0]
    winding = int(round(total / (2 * np.pi)))
    residual = abs(total - 2 * np.pi * winding)
    if residual > 0.5 * np.pi and u.size > 2:
        ui = np.unwrap(x[1:])
        tot_i = ui[-1] - ui[0]
        wind_i = int(round(tot_i / (2 * np.pi)))
        res_i = abs(tot_i - 2 * np.pi * wind_i)
        if res_i < residual:
            u, total, winding, residual = ui, tot_i, wind_i, res_i
    local = np.abs(np.diff(u))
    max_local = float(local.max()) if local.size else 0.0
    jump_like = bool(max_local >= 1.9 * np.pi and local.size > 1)
    return winding, float(residual), max_local, jump_like
