#!/usr/bin/env python3
"""make_win_kpoints (R1.2): generate the wannier90.win 'begin kpoints' block
from the ACTUAL k-points of a finished VASP run (IBZKPT preferred; OUTCAR
fallback).

FORMAT (VASP wiki KPOINTS/IBZKPT, explicit mode - the format VASP writes IBZKPT
in):
  line 1            : comment
  line 2            : number of k-points (nk)
  line 3            : coordinate system ('C/c/K/k' = Cartesian, anything else
                      = fractional/reciprocal). NOTE: VASP's automatic-mesh
                      IBZKPT uses this line for 'Gamma'/'Monkhorst-pack' (a
                      comment) -> k-point lines ALWAYS start at line index 3.
  lines 4+          : kx ky kz weight   (one k-point per line, 4 numbers)

Checks (fail-closed):
  C1 count   : exactly nk k-point lines must be present (indexing fixed: the
               old parser skipped the first k-point by starting at line 5).
  C2 format  : each k-line must have >= 3 coordinates + weight (4 numbers).
  C3 weights : all weights > 0 and sum(weights) == nk within 1e-6 (full-BZ
               mesh, ISYM=-1: every point weight 1).
  C4 order   : first point must be Gamma (0,0,0) for our Gamma-centered mesh.
  C5 coords  : fractional/reciprocal ONLY - reject Cartesian mode lines
               ('C','c','K','k','Cartesian','1') because the win requires
               fractional coordinates. Reject |k_i| > 0.5 + 1e-6.

Usage:
  python3 make_win_kpoints.py IBZKPT [mp_grid "9 9 1"]
  python3 make_win_kpoints.py OUTCAR [mp_grid "9 9 1"]
"""
import sys
from pathlib import Path

MODES_CART = {"c", "k", "cartesian", "1"}
MODES_FRAC = {"r", "f", "reciprocal", "fractional", "0", "gamma", "monkhorst",
              "monkhorst-pack", "mp", "auto"}


def read_ibzkpt(path, nk_expected=None):
    lines = Path(path).read_text(errors="replace").splitlines()
    nk = int(lines[1].split()[0])
    # k-point lines ALWAYS start at index 3 (comment, nk, mode/comment)
    start = 3
    mode = lines[2].strip().lower()
    if mode and mode[0] in MODES_CART:
        raise SystemExit(
            f"FAIL C5: IBZKPT mode '{lines[2].strip()}' is CARTESIAN; the win "
            f"requires fractional/reciprocal coordinates (set KPOINTS to "
            f"'Reciprocal' or an automatic Gamma mesh)")
    klines = [ln.split() for ln in lines[start:start + nk]]
    if len(klines) != nk:
        raise SystemExit(
            f"FAIL C1: found {len(klines)} k-lines, expected {nk} "
            f"(indexing/header mismatch)")
    return _validate(klines, nk, nk_expected)


def read_outcar(path):
    text = Path(path).read_text(errors="replace")
    blocks = []
    cur = None
    for ln in text.splitlines():
        if ln.startswith(" k-points in reciprocal lattice and weights:"):
            cur = []
        elif cur is not None:
            parts = ln.split()
            if len(parts) == 4:
                try:
                    cur.append(([float(x) for x in parts[0:3]], float(parts[3])))
                except ValueError:
                    cur = None
            else:
                cur = None
        if cur is not None and len(cur) > 0 and ln.strip() == "":
            blocks.append(cur)
            cur = None
    if not blocks:
        raise SystemExit("FAIL: no k-point block found in OUTCAR")
    klines = [[f"{k[0]:.12f}", f"{k[1]:.12f}", f"{k[2]:.12f}", f"{w:.12f}"]
              for k, w in blocks[-1]]
    return _validate(klines, None, None)


def _validate(klines, nk, nk_expected):
    n = len(klines)
    if nk is not None and n != nk:
        raise SystemExit(f"FAIL C1: {n} k-lines != declared nk={nk}")
    if nk_expected is not None and n != nk_expected:
        raise SystemExit(f"FAIL C1: {n} k-points != expected {nk_expected}")
    for i, ln in enumerate(klines):
        if len(ln) < 4:
            raise SystemExit(f"FAIL C2: k-line {i}: {ln} (need 3 coords + weight)")
        k = [float(x) for x in ln[0:3]]
        w = float(ln[3])
        if w <= 0:
            raise SystemExit(f"FAIL C3: k-line {i}: non-positive weight {w}")
        if max(abs(x) for x in k) >= 1.0 - 1e-6:
            raise SystemExit(
                f"FAIL C5: k-line {i} outside fractional range [0,1): {k}")
    k0 = [float(x) for x in klines[0][0:3]]
    if not all(abs(x) < 1e-6 for x in k0):
        raise SystemExit(
            f"FAIL C4: first k-point {k0} != Gamma (0,0,0) - "
            f"mesh not Gamma-centered")
    ws = [float(ln[3]) for ln in klines]
    if abs(sum(ws) - n) > 1e-6:
        raise SystemExit(f"FAIL C3: sum(weights)={sum(ws)} != nk={n}")
    return n, [tuple(float(x) for x in ln[0:3]) for ln in klines], ws


def main():
    src = sys.argv[1]
    mp = sys.argv[2] if len(sys.argv) > 2 else None
    try:
        n, ks, ws = read_ibzkpt(src)   # IBZKPT-style header (auto-detect)
    except (SystemExit, IndexError, ValueError) as e:
        if "--force-outcar" in sys.argv:
            n, ks, ws = read_outcar(src)
        else:
            raise SystemExit(
                f"not an IBZKPT-style file ({e}); use --force-outcar for "
                f"OUTCAR parsing")
    print("begin kpoints")
    for k in ks:
        print(" %.12f %.12f %.12f" % k)
    print("end kpoints")
    if mp:
        n1, n2, n3 = (int(x) for x in mp.split())
        exp = n1 * n2 * n3
        print(f"# count={n} (mp_grid {n1}x{n2}x{n3} -> {exp} full-BZ points; "
              f"weights sum={sum(ws):.6f})")


if __name__ == "__main__":
    main()
