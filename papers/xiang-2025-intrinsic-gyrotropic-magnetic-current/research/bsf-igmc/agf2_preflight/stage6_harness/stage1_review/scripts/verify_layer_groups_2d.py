#!/usr/bin/env python3
"""Direct 2D layer-group verification of the AgF2 monolayer and the constructed
B1/B2 bilayer candidates - single-transform setting reconciliation.

Method:
  1. Detect the structure's symmetry operations with spglib (input frame).
  2. Reconcile the detected op set with the pinned toolkit v0.18.0 layer-group
     registry (standard setting of the expected layer group) using ONE global
     origin shift delta (plus, if needed, one in-plane basis permutation S).
     Every operation must map simultaneously with the SAME (S, delta);
     per-operation shifts are rejected as a setting check.
  3. Independently: apply each toolkit operation to the origin-shifted
     structure and require an exact atom-by-atom (same-species, modulo lattice
     in x,y; z compared exactly inside the slab) self-map.

This validates the *2D layer group* of the input structure and is the primary
check; the spglib 3D classification is reported separately.
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

from group_theory_operations import get_crystallographic_layer_group  # noqa: E402
from group_theory_operations import layer_group_operations  # noqa: E402

from gtot_common import (  # noqa: E402
    find_global_setting,
    ops_map_structure,
    read_poscar,
    wrapf,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main() -> None:
    targets = [
        ("monolayer", ROOT / "structures" / "C2DB_2AgF2-1_monolayer" / "POSCAR", 17, "p2_1/b11"),
        ("B1", ROOT / "structures" / "constructed_B1_candidate" / "POSCAR", 12, "pb11"),
        ("B2", ROOT / "structures" / "constructed_B2_candidate" / "POSCAR", 12, "pb11"),
    ]
    report = {}
    for label, path, lg_num, sym in targets:
        lattice, numbers, frac = read_poscar(path)
        ops = layer_group_operations(lg_num)
        detected = []
        import spglib  # noqa: PLC0415

        ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=1e-3)
        detected = [(np.asarray(R, dtype=int), np.asarray(t, dtype=float))
                    for R, t in zip(ds.rotations, ds.translations)]
        recon = find_global_setting(detected, ops)
        delta = recon["delta"]
        mapped, all_ok = ops_map_structure(ops, frac, numbers, delta)
        lg = get_crystallographic_layer_group(lg_num)
        report[label] = {
            "expected_layer_group": f"{lg.international_short} (LG{lg_num})",
            "operation_count_total": len(ops),
            "single_global_transform": {
                "S_basis_permutation": recon["S"],
                "origin_shift_delta_fractional": delta,
                "matched_operations_after_global_shift": recon["matched"],
                "detected_operations": recon["total"],
                "per_operation": recon["detail"],
                "verdict": recon["verdict"],
            },
            "direct_atom_mapping_with_same_delta": {
                "operation_count_matching_structure": sum(1 for r in mapped if r["maps_structure_exactly"]),
                "per_operation": mapped,
                "all_map_exactly": all_ok,
            },
            "verdict": "CONSISTENT_WITH_LAYER_GROUP"
            if recon["verdict"] == "CONSISTENT" and all_ok
            else "SETTING_MISMATCH",
            "notes": "spglib ops reported in the input (POSCAR) frame; the single global "
                     "transform (S, delta) maps the full detected set onto the toolkit "
                     "standard setting simultaneously.",
        }
        print(f"{label}: toolkit {sym} (LG{lg_num}) ops={len(ops)}, "
              f"global-transform reconciled={recon['matched']}/{recon['total']} "
              f"delta={delta}, atom-mapping ok={all_ok} -> {report[label]['verdict']}")
    out = ROOT / "results" / "layer_groups_2d_report.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
