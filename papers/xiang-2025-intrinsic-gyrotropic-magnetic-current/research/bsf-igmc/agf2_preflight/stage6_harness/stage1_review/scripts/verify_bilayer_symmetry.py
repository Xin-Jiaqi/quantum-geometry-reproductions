#!/usr/bin/env python3
"""3D setting reconciliation and spglib classification of the AgF2 monolayer
and the constructed B1/B2 bilayer candidates.

For each structure:
  1. spglib 3D space-group classification of the slab cell (large z vacuum);
  2. reconcile the DETECTED operations (spglib, input frame) against the
     pinned toolkit v0.18.0 layer-group registry (standard setting) with ONE
     global transform (S, delta); every operation must map simultaneously;
  3. cross-reference the direct atom-mapping check of the same toolkit ops
     (see layer_groups_2d_report.json) - both use the same single delta.

Note on 2D vs 3D settings: spglib's *magnetic* symmetry dataset and its
space-group dataset both report rotations in the INPUT frame for these cells;
the toolkit standard setting differs from the input frame by at most the
single origin shift reported here (S = identity for all systems tested).
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

from group_theory_operations import get_crystallographic_layer_group  # noqa: E402
from group_theory_operations import layer_group_operations  # noqa: E402

from gtot_common import find_global_setting, read_poscar  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def analyze(path: Path, expected_lg: int, label: str) -> dict:
    lattice, numbers, frac = read_poscar(path)
    cell = (lattice, frac, numbers)
    report = {"label": label, "poscar": str(path.relative_to(ROOT)), "natoms": len(numbers)}

    ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
    report["spglib"] = {
        "number": int(ds.number),
        "international": ds.international,
        "hall": ds.hall,
        "pointgroup": ds.pointgroup,
        "n_operations": int(len(ds.rotations)),
        "wyckoffs": list(ds.wyckoffs),
        "site_symmetry_symbols": list(ds.site_symmetry_symbols),
        "equivalent_atoms": [int(x) for x in ds.equivalent_atoms],
        "note": "spglib operations reported in the input (POSCAR) frame; "
                "hall symbol labels the conventional space-group type.",
    }
    lg = get_crystallographic_layer_group(expected_lg)
    report["expected_layer_group"] = {
        "number": lg.number,
        "international_short": lg.international_short,
        "international_full": lg.international_full,
        "point_group_hm": lg.point_group_hm,
    }
    ops = layer_group_operations(expected_lg)
    report["expected_layer_group_operation_count"] = len(ops)

    detected = [(np.asarray(R, dtype=int), np.asarray(t, dtype=float))
                for R, t in zip(ds.rotations, ds.translations)]
    recon = find_global_setting(detected, ops)
    report["setting_reconciliation"] = {
        "method": "single global transform (S, delta) applied to ALL detected "
                  "operations simultaneously; per-operation shifts rejected",
        "S_basis_permutation": recon["S"],
        "origin_shift_delta_fractional": recon["delta"],
        "matched_operations_after_global_shift": recon["matched"],
        "detected_operations": recon["total"],
        "per_operation": recon["detail"],
        "verdict": recon["verdict"],
        "cross_reference": "direct atom-mapping with the same delta in "
                           "results/layer_groups_2d_report.json",
    }
    print(f"{label}: spglib SG {ds.international} (#{ds.number}) n_ops={len(detected)}; "
          f"expected LG {lg.international_short} (#{expected_lg}); "
          f"global-transform reconciled ops={recon['matched']}/{recon['total']} "
          f"S={recon['S']} delta={recon['delta']}")
    return report


def main() -> None:
    reports = {}
    reports["monolayer"] = analyze(
        ROOT / "structures" / "C2DB_2AgF2-1_monolayer" / "POSCAR",
        expected_lg=17,
        label="C2DB 2AgF2-1 monolayer",
    )
    reports["B1"] = analyze(
        ROOT / "structures" / "constructed_B1_candidate" / "POSCAR",
        expected_lg=12,  # pb11
        label="constructed B1 candidate",
    )
    reports["B2"] = analyze(
        ROOT / "structures" / "constructed_B2_candidate" / "POSCAR",
        expected_lg=12,  # pb11
        label="constructed B2 candidate",
    )
    out = ROOT / "results" / "bilayer_symmetry_report.json"
    out.write_text(json.dumps(reports, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
