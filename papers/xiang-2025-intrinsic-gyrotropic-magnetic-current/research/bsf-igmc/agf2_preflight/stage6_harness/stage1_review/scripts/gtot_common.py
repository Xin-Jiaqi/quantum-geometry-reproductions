#!/usr/bin/env python3
"""Shared helpers for layer-group / magnetic symmetry verification.

Implements a RIGOROUS single-transform setting reconciliation: one origin
shift delta (fractional, modulo lattice) - optionally combined with one basis
permutation S - must map EVERY operation of the detected set onto the
registry set simultaneously. Per-operation shifts are rejected.

Also provides the direct structure check: with the SAME delta, every toolkit
operation must map the (origin-shifted) structure onto itself atom-by-atom.
"""
from __future__ import annotations

import numpy as np
from pathlib import Path


def wrapf(v):
    """Wrap a fractional coordinate into [0,1) with boundary snapping."""
    w = float(v) % 1.0
    if w > 1.0 - 1e-6:
        return 0.0
    if w < 1e-6:
        return 0.0
    return w


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    lattice = np.array([[float(x) for x in line.split()] for line in lines[2:5]])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coord_type = lines[7].strip().lower()
    natoms = sum(counts)
    coords = np.array([[float(x) for x in lines[8 + i].split()] for i in range(natoms)])
    if coord_type.startswith("car"):
        cart = coords
        frac = np.linalg.solve(lattice.T, cart.T).T
    else:
        frac = coords
    numbers = []
    for sp, c in zip(species, counts):
        z = 47 if sp == "Ag" else 9
        numbers += [z] * c
    return lattice, np.array(numbers), frac


def op_key(R, t):
    return (tuple(map(int, np.asarray(R, dtype=int).ravel())),
            tuple(round(wrapf(x), 6) for x in np.asarray(t, dtype=float)))


def find_global_setting(detected, registry, S_candidates=None, tol=1e-4):
    """Find (S, delta) mapping the detected op set onto the registry op set.

    detected: list of (R, t) ops in the input (POSCAR) frame.
    registry: toolkit SeitzOp list (standard setting of the expected layer
              group).
    S_candidates: list of 3x3 basis permutation matrices to try (identity plus
              axis relabelings of the in-plane axes).

    Returns dict with S, delta, matched, total, verdict and a per-op mapping
    table. delta is ONE global fractional origin shift (modulo lattice and
    modulo the kernel of (S R S^-1 - I)); per-op shifts are not used.
    """
    if S_candidates is None:
        S_candidates = [np.eye(3, dtype=int)]
    reg_by_key = {op_key(op.rotation, op.translation) for op in registry}
    best = None
    for S in S_candidates:
        S = np.asarray(S, dtype=int)
        # rotate detected ops into the registry frame under S
        trans = []
        for R, t in detected:
            Rs = S @ R @ S.T
            ts = S @ t
            trans.append((Rs, ts))
        # rotation sets must match exactly
        det_rot = {tuple(map(int, Rs.ravel())) for Rs, _ in trans}
        reg_rot = {tuple(map(int, np.asarray(op.rotation, dtype=int).ravel())) for op in registry}
        if det_rot != reg_rot:
            continue
        # build equations (S R S^-1 - I) delta = t_reg - S t  (mod 1) over ops with R != I
        A = []
        b = []
        pairs = []
        for (Rs, ts), (R, t) in zip(trans, detected):
            # match registry op with same rotation
            regs = [op for op in registry
                    if tuple(map(int, np.asarray(op.rotation, dtype=int).ravel()))
                    == tuple(map(int, Rs.ravel()))]
            if not regs:
                return None
            reg_op = regs[0]
            pairs.append(((Rs, ts), (np.asarray(reg_op.rotation, dtype=int),
                                     np.asarray(reg_op.translation, dtype=float))))
            M = Rs - np.eye(3)
            if np.any(M != 0):
                A.append(M)
                b.append((reg_op.translation - ts) % 1.0)
        if not A:
            delta = np.zeros(3)
        else:
            A = np.vstack(A)
            b = np.concatenate(b)
            delta = np.linalg.pinv(A) @ b
            delta = np.array([wrapf(x) for x in delta])
        # verify ALL ops with the SAME delta
        matched = 0
        detail = []
        for (Rs, ts), (Rreg, treg) in pairs:
            t2 = ts + (Rs - np.eye(3)) @ delta
            ok = op_key(Rs, t2) == op_key(Rreg, treg)
            matched += int(ok)
            detail.append({"rotation": tuple(map(int, Rs.ravel())),
                           "detected_translation": [round(x, 6) for x in ts],
                           "mapped_translation": [round(wrapf(x), 6) for x in t2],
                           "registry_translation": [round(wrapf(x), 6) for x in treg],
                           "match": bool(ok)})
        if best is None or matched > best["matched"]:
            best = {"S": S.tolist(), "delta": [round(x, 6) for x in delta],
                    "matched": matched, "total": len(detected), "detail": detail}
        if matched == len(detected):
            break
    if best is None:
        return {"verdict": "SETTING_MISMATCH", "S": None, "delta": None,
                "matched": 0, "total": len(detected), "detail": []}
    best["verdict"] = ("CONSISTENT" if best["matched"] == best["total"]
                       else "SETTING_MISMATCH")
    return best


def ops_map_structure(ops, frac, numbers, delta=None):
    """Check each Seitz op maps the structure onto itself in ONE fixed frame.

    The structure is first shifted by the single global delta (origin
    transformation into the toolkit standard setting); then every op must be an
    exact symmetry (atom-by-atom, same species, modulo the 3D periodic cell).
    z is compared modulo the cell height: the slab lives in a periodic box with
    large vacuum, and layer-group ops that flip z (inversion, screw axes) map
    the slab onto its periodic image; this matches spglib's 3D classification
    of the same cell. The 2D (layer-group) statement follows from the full op
    set acting on the layer. All coordinates are wrapped with boundary
    snapping; a single delta is used for all ops (per-op shifts rejected)."""
    delta = np.zeros(3) if delta is None else np.array(delta)
    natoms = len(numbers)
    results = []
    all_ok = True
    for op in ops:
        R = np.asarray(op.rotation, dtype=float)
        t = np.asarray(op.translation, dtype=float)
        used = {}
        for i in range(natoms):
            v = frac[i] + delta
            key = (round(wrapf(v[0]), 6), round(wrapf(v[1]), 6), round(wrapf(v[2]), 6),
                   int(numbers[i]))
            used[key] = used.get(key, 0) + 1
        ok = True
        for i in range(natoms):
            v = frac[i] + delta
            vp = R @ v + t
            key = (round(wrapf(vp[0]), 6), round(wrapf(vp[1]), 6), round(wrapf(vp[2]), 6),
                   int(numbers[i]))
            if used.get(key, 0) == 0:
                ok = False
                break
            used[key] -= 1
        all_ok = all_ok and ok
        results.append({"rotation": np.asarray(R, dtype=int).ravel().tolist(),
                        "translation": [round(x, 6) for x in t],
                        "maps_structure_exactly": bool(ok)})
    return results, all_ok
