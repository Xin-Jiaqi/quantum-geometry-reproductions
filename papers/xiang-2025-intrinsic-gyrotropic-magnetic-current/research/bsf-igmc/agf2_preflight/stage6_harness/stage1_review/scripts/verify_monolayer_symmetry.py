#!/usr/bin/env python3
"""Verify the C2DB AgF2 monolayer structure:
  1. spglib 3D space-group classification of the slab cell;
  2. layer-group operations p2_1/b11 (LG17) from the pinned
     group-theory-operations-toolkit (v0.18.0) compared with the detected
     operations, in the monolayer's own coordinate setting.

This check is independent of any DFT output; it validates the *input*
structure only.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.environ.get("GTOT_TOOLKIT_SRC", "/tmp/gtot/toolkit-v0.18.0/src"))
# (toolkit source path is configurable; the pinned clone is at /tmp/gtot/toolkit-v0.18.0,
# tag v0.18.0 = commit 40a68e9ee2235a765e768f07f0f11fb99729c5c2, repo Xin-Jiaqi/group-theory-operations-toolkit)
import spglib  # noqa: E402

from group_theory_operations import layer_group_operations  # noqa: E402
from group_theory_operations import get_crystallographic_layer_group  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
POSCAR = ROOT / "structures" / "C2DB_2AgF2-1_monolayer" / "POSCAR"


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    lattice = np.array([[float(x) for x in line.split()] for line in lines[2:5]])
    lattice = lattice * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coord_type = lines[7].strip().lower()
    natoms = sum(counts)
    coords = []
    for i in range(8, 8 + natoms):
        coords.append([float(x) for x in lines[i].split()])
    coords = np.array(coords)
    if coord_type.startswith("car") or coord_type.startswith("cart"):
        cart = coords
        frac = np.linalg.solve(lattice.T, cart.T).T
    else:
        frac = coords
        cart = (frac @ lattice)
    numbers = []
    for sp, c in zip(species, counts):
        z = 47 if sp == "Ag" else 9
        numbers += [z] * c
    return lattice, np.array(numbers), cart, frac


def main() -> None:
    lattice, numbers, cart, frac = read_poscar(POSCAR)
    cell = (lattice, frac, numbers)

    report = {"structure": "C2DB 2AgF2-1 monolayer",
              "poscar": str(POSCAR.relative_to(ROOT))}

    for symprec in (1e-4, 1e-3):
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        report[f"spglib_symprec_{symprec:g}"] = {
            "number": int(ds.number),
            "international": ds.international,
            "hall": ds.hall,
            "pointgroup": ds.pointgroup,
            "n_operations": int(len(ds.rotations)),
        }
        print(f"symprec={symprec}: SG {ds.international} (#{ds.number}), "
              f"n_ops={len(ds.rotations)}, pg={ds.pointgroup}")

    # layer-group check: operations of LG17 p2_1/b11 in toolkit registry
    lg = get_crystallographic_layer_group(17)
    report["toolkit_layer_group"] = {
        "number": lg.number,
        "international_short": lg.international_short,
        "international_full": lg.international_full,
        "crystal_system": lg.crystal_system,
        "point_group_hm": lg.point_group_hm,
    }
    print("toolkit LG17:", lg.international_short, lg.international_full)

    ops = layer_group_operations(17)
    ops_list = []
    for op in ops:
        ops_list.append({"rotation": op.rotation.tolist() if hasattr(op.rotation, "tolist") else list(op.rotation),
                         "translation": list(op.translation) if hasattr(op.translation, "__iter__") else op.translation})
    report["toolkit_lg17_operation_count"] = len(ops_list)

    # Compare with spglib-detected operations: ONE global origin shift must map
    # the whole detected op set onto the toolkit LG17 standard setting
    # simultaneously (per-operation shifts are rejected).
    from gtot_common import find_global_setting, ops_map_structure  # noqa: PLC0415

    ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
    detected = [(np.asarray(R, dtype=int), np.asarray(t, dtype=float))
                for R, t in zip(ds.rotations, ds.translations)]
    recon = find_global_setting(detected, ops)
    mapped, all_ok = ops_map_structure(ops, frac, numbers, recon["delta"])
    report["detected_vs_lg17_matches"] = {
        "detected": recon["total"],
        "matched_to_lg17_after_single_global_transform": recon["matched"],
        "S_basis_permutation": recon["S"],
        "origin_shift_delta_fractional": recon["delta"],
        "per_operation": recon["detail"],
        "direct_atom_mapping_all_ops_exact": all_ok,
        "verdict": recon["verdict"],
    }
    print(f"detected ops={recon['total']}, matched to toolkit LG17 ops="
          f"{recon['matched']} (single global delta={recon['delta']}), "
          f"atom-mapping ok={all_ok}")

    out = ROOT / "results" / "monolayer_symmetry_report.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
