#!/usr/bin/env python3
"""make_win_kpoints: generate the wannier90.win 'begin kpoints' block from the
ACTUAL k-points of a finished VASP run (IBZKPT preferred; OUTCAR fallback).

Why: the VASP-W90 interface requires the win k-point list to match, point for
point, the k-points VASP actually computed (order + coordinates). With ISYM=-1
(no symmetry reduction, as set in the soc_w90 INCARs) a Gamma-centered 9x9x1
mesh runs the full 81 k-points. Never hand-edit the kpoints block.

Usage:
  python3 make_win_kpoints.py IBZKPT [mp_grid "9 9 1"] > kpoints_block.txt
  python3 make_win_kpoints.py OUTCAR [mp_grid "9 9 1"] > kpoints_block.txt
Then replace the placeholder line in wannier90.win with the generated block:
  # begin kpoints: GENERATED POST-SCF ...
"""
import sys
from pathlib import Path


def read_ibzkpt(path):
    lines = Path(path).read_text().splitlines()
    # IBZKPT: comment, nk, comment, cartesian flag(0), then nk lines kx ky kz w
    nk = int(lines[1].split()[0])
    ks, ws = [], []
    for ln in lines[4:4 + nk]:
        parts = ln.split()
        ks.append(tuple(float(x) for x in parts[0:3]))
        ws.append(float(parts[3]))
    return ks, ws


def read_outcar(path):
    text = Path(path).read_text()
    # fallback: parse the "k-points in reciprocal lattice and weights" blocks;
    # take the LAST block (final SCF/interface step).
    blocks = []
    cur = None
    for ln in text.splitlines():
        if ln.startswith(" k-points in reciprocal lattice and weights:"):
            cur = []
        elif cur is not None:
            parts = ln.split()
            if len(parts) == 4:
                try:
                    cur.append((tuple(float(x) for x in parts[0:3]), float(parts[3])))
                except ValueError:
                    cur = None
            else:
                cur = None
        if cur is not None and len(cur) > 0 and ln.strip() == "":
            blocks.append(cur)
            cur = None
    if not blocks:
        raise SystemExit("no k-point block found in OUTCAR")
    return blocks[-1]


def main():
    src = sys.argv[1]
    mp = sys.argv[2] if len(sys.argv) > 2 else None
    if Path(src).name.upper() == "IBZKPT":
        ks, ws = read_ibzkpt(src)
    else:
        pairs = read_outcar(src)
        ks, ws = zip(*pairs)
    print(f"begin kpoints")
    for k in ks:
        print(" %.12f %.12f %.12f" % k)
    print("end kpoints")
    if mp:
        n1, n2, n3 = (int(x) for x in mp.split())
        print(f"# count={len(ks)} (mp_grid {n1}x{n2}x{n3} -> {n1*n2*n3} full-BZ points)")


if __name__ == "__main__":
    main()
