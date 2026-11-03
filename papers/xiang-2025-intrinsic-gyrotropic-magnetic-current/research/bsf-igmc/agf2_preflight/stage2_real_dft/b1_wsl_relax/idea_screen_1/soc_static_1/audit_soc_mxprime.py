#!/usr/bin/env python3
"""Independent Ag4F8 SOC spin/orbital moment m_x' compatibility audit.

Usage: python audit_soc_mxprime.py MOMENTS.csv
Read-only, Python stdlib. Pair mapping was independently verified against the
actual POSCAR glide in the earlier screening structural audit; the current
script does not re-derive geometry. Prints JSON, no file modification.

m_x' antiunitary transforms time-odd axial M=(Mx,My,Mz) to (-Mx,My,Mz).
This tests *site-projected output precision* only; it cannot certify full
magnetic space group, density-grid symmetries, easy axis or orbital IGMC.
"""
import csv
import json
import math
import sys

PAIRS = [(1,2),(3,4),(5,8),(6,7),(9,12),(10,11)]
AXES = ["x", "y", "z"]
T = [-1.0,1.0,1.0]

def compute(path):
    with open(path, newline="", encoding="utf-8") as src:
        rows = {int(row["atom_1based"]): row for row in csv.DictReader(src)}
    if sorted(rows) != list(range(1,13)):
        raise ValueError("expected 12 Ag4F8 rows in fixed order")
    if any(rows[i]["species"] != ("Ag" if i <=4 else "F") for i in range(1,13)):
        raise ValueError("unexpected Ag4F8 species ordering")
    checks = []
    for channel in ("S", "L"):
        for i,j in PAIRS:
            ai = [float(rows[i][channel+a+"_mub"]) for a in AXES]
            aj = [float(rows[j][channel+a+"_mub"]) for a in AXES]
            expected = [T[k]*ai[k] for k in range(3)]
            err = [aj[k]-expected[k] for k in range(3)]
            checks.append({"channel":channel,"pair":[i,j],
                           "expected_destination_mub":expected,
                           "observed_destination_mub":aj,
                           "residual_mub":err,
                           "max_abs_residual_mub":max(map(abs,err))})
    sums = {}
    for channel in ("S","L"):
        sums[channel] = [sum(float(rows[i][channel+a+"_mub"]) for i in range(1,13)) for a in AXES]
    return {
        "evidence": "PR #38 committed site-projected spin and orbital CSV only",
        "operation": "antiunitary m_x' maps (Mx,My,Mz) -> (-Mx,+My,+Mz)",
        "source_structure_site_pairs_1based": PAIRS,
        "max_spin_abs_residual_mub": max(c["max_abs_residual_mub"] for c in checks if c["channel"]=="S"),
        "max_orbital_abs_residual_mub": max(c["max_abs_residual_mub"] for c in checks if c["channel"]=="L"),
        "output_rounding_tolerance_mub": 0.00101,
        "passes_reported_precision": all(c["max_abs_residual_mub"] <= 0.00101 for c in checks),
        "sums_from_rounded_projected_moments_mub":sums,
        "caution": "Projected site sums are not automatically total bulk spin/orbital moments; not full MPG certification.",
        "pairs":checks,
    }

def main():
    if len(sys.argv) !=2:
        raise SystemExit("Usage: python audit_soc_mxprime.py MOMENTS.csv")
    print(json.dumps(compute(sys.argv[1]),indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
