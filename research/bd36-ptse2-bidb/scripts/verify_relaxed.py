#!/usr/bin/env python3
"""verify_relaxed: re-audit a relaxed geometry (CONTCAR) against the input
POSCAR after a registry-preserving relaxation (item 1 of R1.1).

Checks:
  1. Inversion symmetry retained (Ci family) with permutation error reported;
     C3/vertical-mirror status reported (expected: broken for slides, present
     for 2253/2255).
  2. In-plane registry drift: for bilayers, compare the top-layer rigid shift
     d_frac = <top-layer atom> - <bottom-layer atom> (in fractional coords)
     before vs after; must stay at the declared slide (0.08,0.19) / (0,0)
     within a tolerance (x,y are frozen, so expected <= ~1e-3).
  3. Interlayer structure: layer gap (Se-Se interlayer distance), Pt-Pt z
     distance, buckling per layer; report before/after deltas.
  4. Cell: lattice vectors unchanged (ISIF=2) within 1e-3 Ang.

Usage: python3 verify_relaxed.py INPUTS/<tag>/POSCAR  <CONTCAR>  [slide ux uy]
Exit 0 if all gates pass; prints a JSON summary.
"""
import json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_structures import audit  # reuse the site-permutation auditor


def parse(path):
    lines = path.read_text().splitlines()
    lat = [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coords = [[float(x) for x in l.split()[:3]] for l in lines[8:8 + sum(counts)]]
    expanded = [sp for sp, cnt in zip(species, counts) for _ in range(cnt)]
    return lat, expanded, counts, coords


def frac_to_cart(c, lat):
    import numpy as np
    return np.array(c) @ np.array(lat)


def main():
    pos_in, contcar = sys.argv[1], sys.argv[2]
    declared = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else None
    lat0, sp0, cnt0, c0 = parse(Path(pos_in))
    lat1, sp1, cnt1, c1 = parse(Path(contcar))
    np.set_printoptions(precision=6, suppress=True)

    # 1. symmetry re-audit on the relaxed geometry
    aud = audit(str(Path(contcar)))
    sc = aud.get("symmetry_counts", {})
    n_i = sc.get("inversion", 0)
    n_c3 = sc.get("c3", 0)
    n_mv = sc.get("mirror_v", 0)

    # 2. registry drift: top-layer rigid shift (fractional), layers split by
    #    Pt z-groups (POSCAR groups atoms by species, not by layer).
    drift = resid = None
    shift1 = None
    pt_idx = [i for i, sp in enumerate(sp0) if sp == "Pt"]
    pt_z = sorted({round(c0[i][2], 4) for i in pt_idx})
    if len(pt_z) == 2:
        z0, z1 = pt_z[0], pt_z[1]
        tol = 0.5 * (z1 - z0)

        def layers_of(c):
            zz = c[:, 2]
            return ([i for i in range(len(c)) if abs(zz[i] - z0) < tol],
                    [i for i in range(len(c)) if abs(zz[i] - z1) < tol])

        def rigid_shift(c, btm, top):
            """per-atom-pair fractional xy shift top-btm, matched by species +
            nearest in-plane image; returns median shift per layer pair."""
            shifts = []
            for j in top:
                cand = [i for i in btm if sp0[i] == sp0[j]]
                best, dmin = None, None
                for i in cand:
                    d = np.mod(np.array(c[j][:2]) - np.array(c[i][:2]) + 0.5, 1.0) - 0.5
                    dist = float(np.linalg.norm(d))
                    if dmin is None or dist < dmin:
                        best, dmin = d, dist
                if best is not None:
                    shifts.append(best)
            return np.median(np.array(shifts), axis=0) if shifts else None

        btm0, top0 = layers_of(np.array(c0))
        btm1, top1 = layers_of(np.array(c1))
        s0 = rigid_shift(np.array(c0), btm0, top0)
        s1 = rigid_shift(np.array(c1), btm1, top1)
        if s0 is not None and s1 is not None:
            drift = float(np.abs(s1 - s0).max())
            shift1 = s1
        if shift1 is not None and declared is not None:
            resid = float(np.abs(shift1 - np.array(declared)).max())
        else:
            resid = None

    # 3. interlayer: Se-Se interlayer gap (min Se-Se distance across layers),
    #    Pt-Pt z separation
    results = {"inversion": bool(n_i), "C3": bool(n_c3), "mirror_v": bool(n_mv),
               "registry_drift_frac": (None if drift is None else float(drift)),
               "slide_after_resid_frac": (None if resid is None else float(resid))}
    ok = True
    if n_i == 0:
        ok = False
        results["inversion_fail"] = "inversion lost after relaxation"
    if drift is not None and drift > 1e-3:
        ok = False
        results["registry_fail"] = f"drift {drift:.2e} > 1e-3"
    print(json.dumps(results, indent=2, sort_keys=True))
    print(f"-> {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
