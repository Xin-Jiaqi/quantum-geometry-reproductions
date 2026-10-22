#!/usr/bin/env python3
"""Read-only Ag4F8 screen-stage symmetry and tensor sign audit (stdlib).

Usage: python check_magnetic_symmetry_screen.py POSCAR.screen_last_complete
No VASP, POTCAR, external package, network or write operation is required.
Uses the prior, explicit collinear spin-axis-x assumption. Magnetic unitary
and antiunitary connector classification is conditional on that assumption.
"""
import math
import sys
from pathlib import Path

SIGNS = (-1, 1, 1, -1)
AXES = "xyz"
R_MX = (-1, 1, 1)
R_MZ = (1, 1, -1)

def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale = float(lines[1])
    if scale <= 0:
        raise ValueError("positive POSCAR scale required")
    basis = [[scale * float(t) for t in lines[i].split()[:3]] for i in range(2, 5)]
    if any(abs(basis[i][j]) > 1e-8 for i in range(3) for j in range(3) if i != j):
        raise ValueError("expected frozen orthorhombic cell")
    cell = [basis[i][i] for i in range(3)]
    if lines[5].split() != ["Ag", "F"] or lines[6].split() != ["4", "8"]:
        raise ValueError("expected Ag4F8 order")
    k = 7
    if lines[k].strip().lower().startswith("s"):
        k += 1
    direct = lines[k].strip().lower().startswith("d")
    if not (direct or lines[k].strip().lower().startswith(("c", "k"))):
        raise ValueError("expected direct/cartesian positions")
    pos = []
    for line in lines[k+1:k+13]:
        p = [float(x) for x in line.split()[:3]]
        pos.append([p[j] * (cell[j] if direct else scale) for j in range(3)])
    if len(pos) != 12:
        raise ValueError("need 12 positions")
    return cell, pos

def nearest_image(v, L):
    return (v + L/2) % L - L/2

def distance(a, b, cell):
    dr = [nearest_image(a[j] - b[j], cell[j]) if j<2 else a[j]-b[j] for j in range(3)]
    return math.sqrt(sum(x*x for x in dr))

def find_mx_glide(cell, pos):
    a, b, c = cell
    mapping = []
    used = set()
    for i, p in enumerate(pos):
        target = [-p[0] + a/2, p[1] + b/2, p[2]]
        pool = range(0, 4) if i < 4 else range(4, 12)
        d, j = min((distance(target, pos[j], cell), j) for j in pool)
        if j in used:
            raise AssertionError("site correspondence not bijective")
        used.add(j)
        mapping.append((i+1, j+1, d))
        if i < 4 and SIGNS[j] != -SIGNS[i]:
            raise AssertionError("m_x' does not map seed spins antiunitarily")
    return mapping

def transform_sign(R, anti, i, j):
    det = math.prod(R)
    return ((-1 if anti else 1) * det * R[i] * R[j])

def tensor_audit():
    allowed = []
    forbidden = []
    rows = []
    for i in range(3):
        for j in range(3):
            name = "alpha_" + AXES[i] + AXES[j]
            same = transform_sign(R_MX, True, i, j)
            sign1 = transform_sign(R_MZ, False, i, j)
            sign2 = transform_sign(R_MZ, True, i, j)
            rows.append((name, "allowed" if same==1 else "forbidden",
                         "even" if sign1==1 else "odd",
                         "even" if sign2==1 else "odd"))
            (allowed if same==1 else forbidden).append(name)
    assert set(allowed)=={"alpha_xx","alpha_yy","alpha_zz","alpha_yz","alpha_zy"}
    assert rows[0]==("alpha_xx","allowed","odd","even")
    assert rows[4]==("alpha_yy","allowed","odd","even")
    return rows

def main():
    if len(sys.argv)!=2:
        raise SystemExit("usage: python check_magnetic_symmetry_screen.py POSCAR")
    cell,pos=read_poscar(sys.argv[1])
    mapping=find_mx_glide(cell,pos)
    maxerr=max(d for _,_,d in mapping)
    print("AgF2 B1 frozen screening snapshot: 12/12 m_x glide site mapping")
    print("operation: (x,y,z) -> (-x + a/2, y + b/2, z)")
    print("max site residual: %.8f A" % maxerr)
    for i,j,d in mapping:
        print("  %02d -> %02d   mismatch %.7g A" % (i,j,d))
    print("Ag x-axis signs assumed:",SIGNS)
    print("Ag magnetic m_x' sign check: PASS")
    print("alpha=P_polar/B_axial, m_x' internal constraint and domain parity:")
    for name,allowed,typeI,typeII in tensor_audit():
        print("  %-10s %-10s  B2/B1 Type-I m_z: %-4s | Type-II m_z': %s" %
              (name,allowed,typeI,typeII))
    print("Pz(B2)=-Pz(B1) for either m_z or m_z' *if a physical domain mapping applies*")
    print("AHE_z parity: m_z even, m_z' odd (allowedness != nonzero).")
    print("All tensor rules are symmetry-only. B2 not relaxed or calculated.")
    if maxerr>1e-2:
        raise AssertionError("screen structure no longer matches m_x glide within 0.01 A")

if __name__ == "__main__":
    main()
