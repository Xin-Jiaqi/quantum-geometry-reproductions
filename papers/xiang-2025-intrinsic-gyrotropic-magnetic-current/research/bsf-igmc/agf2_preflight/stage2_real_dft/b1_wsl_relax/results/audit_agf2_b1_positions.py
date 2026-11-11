#!/usr/bin/env python3
"""Independent read-only Ag4F8 B1 coordinate audit, no VASP execution.

Usage:
  python audit_agf2_b1_positions.py path/to/initial_POSCAR path/to/POSCAR.relaxed
The atom mapping is the frozen C2DB-derived B1 construction in this repository.
Output is JSON to stdout; no input file is changed.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

TOP_AG = (0, 1)
BOTTOM_AG = (2, 3)
TOP_F = (4, 5, 6, 7)
BOTTOM_F = (8, 9, 10, 11)
TOP_INNER_F = (4, 7)
BOTTOM_INNER_F = (8, 11)

def read_poscar(filename):
    data = Path(filename).read_bytes()
    lines = data.decode("utf-8").splitlines()
    scale = float(lines[1])
    if scale <= 0:
        raise ValueError("Expected positive POSCAR scale")
    vectors = [[scale * float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]
    if any(abs(vectors[i][j]) > 1e-8 for i in range(3) for j in range(3) if i != j):
        raise ValueError("Expected orthorhombic B1 lattice")
    cell = [vectors[i][i] for i in range(3)]
    if min(cell) <= 0 or lines[5].split() != ["Ag", "F"] or lines[6].split() != ["4", "8"]:
        raise ValueError("Expected Ag4F8 with positive a,b,c")
    j = 7
    if lines[j].strip().lower().startswith("s"):
        j += 1
    mode = lines[j].strip().lower()
    if not (mode.startswith("c") or mode.startswith("k") or mode.startswith("d")):
        raise ValueError("Unknown POSCAR position mode")
    xyz = []
    for line in lines[j + 1:j + 13]:
        nums = [float(x) for x in line.split()[:3]]
        if len(nums) != 3:
            raise ValueError("Expected 12 complete atomic positions")
        if mode.startswith("d"):
            nums = [nums[i] * cell[i] for i in range(3)]
        else:
            nums = [x * scale for x in nums]
        xyz.append(nums)
    if len(xyz) != 12:
        raise ValueError("Expected 12 atomic positions")
    return {"sha256": hashlib.sha256(data).hexdigest(), "cell": cell, "xyz": xyz}

def mean_z(xyz, inds):
    return sum(xyz[i][2] for i in inds) / len(inds)

def minimum_image_difference(a, b, lattice):
    return (a - b + lattice / 2) % lattice - lattice / 2

def metrics(data):
    x, cell = data["xyz"], data["cell"]
    best = None
    for i in TOP_F:
        for j in BOTTOM_F:
            dx = minimum_image_difference(x[i][0], x[j][0], cell[0])
            dy = minimum_image_difference(x[i][1], x[j][1], cell[1])
            dz = x[i][2] - x[j][2]  # no z-image across the vacuum
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            if best is None or dist < best["distance_A"]:
                best = {"distance_A": dist, "top_ion": i + 1,
                        "bottom_ion": j + 1, "delta_xyz_A": [dx, dy, dz]}
    return {
        "top_Ag_z_A": mean_z(x, TOP_AG),
        "bottom_Ag_z_A": mean_z(x, BOTTOM_AG),
        "Ag_layer_separation_z_A": mean_z(x, TOP_AG) - mean_z(x, BOTTOM_AG),
        "top_inner_F_z_A": mean_z(x, TOP_INNER_F),
        "bottom_inner_F_z_A": mean_z(x, BOTTOM_INNER_F),
        "vertical_inner_F_plane_gap_A": mean_z(x, TOP_INNER_F) - mean_z(x, BOTTOM_INNER_F),
        "minimum_3D_interlayer_FF": best,
    }

def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python audit_agf2_b1_positions.py START_POSCAR STEP15_POSCAR")
    old, new = map(read_poscar, sys.argv[1:])
    if any(abs(x - y) > 1e-6 for x, y in zip(old["cell"], new["cell"])):
        raise ValueError("Fixed-cell input must match")
    cell = old["cell"]
    displacements = [
        [minimum_image_difference(new["xyz"][i][j], old["xyz"][i][j], cell[j])
         if j < 2 else new["xyz"][i][j] - old["xyz"][i][j]
         for j in range(3)]
        for i in range(12)
    ]
    top_dy = sum(displacements[i][1] for i in TOP_AG) / 2
    bottom_dy = sum(displacements[i][1] for i in BOTTOM_AG) / 2
    print(json.dumps({
        "scope": "frozen-coordinate-only; no VASP force/relaxation validation",
        "initial_source_sha256": old["sha256"],
        "step15_source_sha256": new["sha256"],
        "initial": metrics(old),
        "step15": metrics(new),
        "Ag_top_y_translation_A": top_dy,
        "Ag_bottom_y_translation_A": bottom_dy,
        "relative_Ag_layer_y_slide_A": top_dy - bottom_dy,
        "relative_Ag_layer_y_slide_fraction_b": (top_dy - bottom_dy) / cell[1],
        "ion_displacements_minimum_xy_image_A": displacements,
    }, indent=2))

if __name__ == "__main__":
    main()
