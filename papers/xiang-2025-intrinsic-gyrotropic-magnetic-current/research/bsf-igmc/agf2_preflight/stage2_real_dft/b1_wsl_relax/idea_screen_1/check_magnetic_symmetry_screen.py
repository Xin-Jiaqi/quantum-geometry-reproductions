#!/usr/bin/env python3
"""Read-only Ag4F8 structural glide, conditional magnetic pairing and IGMC tensor audit.

Usage: python check_magnetic_symmetry_screen.py POSCAR.screen_last_complete
Python standard library only; no VASP, POTCAR, network or file writes.

Assumptions:
- orthorhombic frozen cell and Ag4 F8 POSCAR order;
- collinear site signs [-,+,+,-] are seeded along physical x;
- m_x glide is {m_x|(1/2,1/2,0)} in the supplied POSCAR origin;
- channel-resolved C/D symmetry protocol from the repository's tensor_symmetry.py.
This is a symmetry audit, not a response calculation or proof of the SOC easy axis.
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
    dr = [nearest_image(a[j] - b[j], cell[j]) if j < 2 else a[j] - b[j] for j in range(3)]
    return math.sqrt(sum(x*x for x in dr))

def find_mx_glide(cell, pos):
    a, b, _ = cell
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
            raise AssertionError("m_x' does not map the collinear seed spins")
    return mapping

def transform_sign(diagonal_R, antiunitary, i, j, channel):
    """Sign for a polar-current / axial-field rank-2 tensor."""
    det = math.prod(diagonal_R)
    eta_T = -1 if (antiunitary and channel == "D") else 1
    return eta_T * det * diagonal_R[i] * diagonal_R[j]

def allowed_components(R, antiunitary, channel):
    return [AXES[i] + AXES[j] for i in range(3) for j in range(3)
            if transform_sign(R, antiunitary, i, j, channel) == 1]

def connector_parity(R, antiunitary, channel, i, j):
    return "same" if transform_sign(R, antiunitary, i, j, channel) == 1 else "reverse"

def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python check_magnetic_symmetry_screen.py POSCAR")
    cell, pos = read_poscar(sys.argv[1])
    mapping = find_mx_glide(cell, pos)
    maxerr = max(d for _, _, d in mapping)
    print("AgF2 B1 screening snapshot: 12/12 m_x glide site mapping")
    print("operation: (x,y,z) -> (-x + a/2, y + b/2, z)")
    print("max site residual: %.8f A" % maxerr)
    for i, j, d in mapping:
        print("  %02d -> %02d   mismatch %.7g A" % (i, j, d))
    print("Ag seed signs along assumed physical x:", SIGNS)
    print("Ag magnetic m_x' pairing check: PASS")
    print("Conditional internal MPG: m_x' = {E,m_x'}; SOC easy axis is NOT established")
    print("Allowed single-state IGMC components under channel-resolved protocol:")
    print("  C (Zeeman-Berry-curvature / conduction):", ", ".join(allowed_components(R_MX, True, "C")))
    print("  D (Zeeman-metric / displacement):", ", ".join(allowed_components(R_MX, True, "D")))
    print("  Static alpha=P/B has the D-like time parity, not the full IGMC tensor.")
    print("B1 -> B2 connector parity (source-defined):")
    print("  component   Type-I m_z C   Type-II m_z' C   Type-I m_z D   Type-II m_z' D")
    for name in ("xy", "yx", "xz", "zx", "xx", "yy", "zz", "yz", "zy"):
        i, j = AXES.index(name[0]), AXES.index(name[1])
        c1 = connector_parity(R_MZ, False, "C", i, j)
        c2 = connector_parity(R_MZ, True, "C", i, j)
        d1 = connector_parity(R_MZ, False, "D", i, j)
        d2 = connector_parity(R_MZ, True, "D", i, j)
        print("  %-10s %-15s %-16s %-15s %s" % (name, c1, c2, d1, d2))
    print("Pz reverses for both structural partners. AHE sigma_xy: Type-I same; Type-II reverses.")
    print("These are symmetry-only predictions; actual SOC MPG and response amplitudes are uncomputed.")
    if maxerr > 1e-2:
        raise AssertionError("screen structure no longer matches m_x glide within 0.01 A")

if __name__ == "__main__":
    main()
