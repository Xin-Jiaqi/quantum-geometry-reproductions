#!/usr/bin/env python3
"""verify_relaxed (R1.2): re-audit a relaxed geometry (CONTCAR, with or without
the Selective-dynamics header) against the input POSCAR.

Checks:
  1. Symmetry re-audit (inversion / C3 / vertical mirrors) via audit_structures.
  2. Per-atom displacement decomposition in CARTESIAN coordinates:
       dr_cart = cart(relaxed) - cart(input)   (minimal in-plane image)
       in-plane |dr_xy|  : must be ~0 for registry-preserving bilayers
                           (x,y FIXED via 'F F T'); reported for all atoms.
       out-of-plane dr_z : the actual relaxation; reported (no hard gate).
     Monolayer (T T T): no in-plane gate.
  3. Registry: top-layer rigid in-plane shift vs input; slide residual vs the
     DECLARED slide (0.08,0.19) when given.
  4. Interlayer: layer gap / Pt-Pt z (reported via the audit).

Selftest (simulated relaxations): python3 verify_relaxed.py --selftest
  - fixture A: z-only displacement (+0.02 A on top layer, -0.01 A on bottom):
    in-plane ~0, dr_z = expected, inversion retained -> PASS
  - fixture B: in-plane drift 0.05 A on one atom (simulates a broken registry):
    must FAIL (in-plane gate)
Usage: python3 verify_relaxed.py INPUT_POSCAR CONTCAR [slide_ux slide_uy]
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_structures import audit


def parse(path):
    """Handle POSCAR/CONTCAR with or without 'Selective dynamics' header.
    Returns (lat, species_expanded, counts, coords_frac, flags_or_None)."""
    lines = path.read_text().splitlines()
    lat = [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    sd = lines[7].strip().lower().startswith("selective")
    start = 9 if sd else 8
    coords = []
    flags = []
    for l in lines[start:start + n]:
        p = l.split()
        coords.append([float(x) for x in p[0:3]])
        if sd and len(p) >= 6:
            flags.append(p[3:6])
    if len(coords) != n:
        raise SystemExit(f"FAIL: {len(coords)} atom lines, expected {n}")
    expanded = [sp for sp, cnt in zip(species, counts) for _ in range(cnt)]
    return (np.array(lat), expanded, counts, np.array(coords),
            flags if sd and len(flags) == n else None)


def frac_to_cart(f, lat):
    return np.asarray(f) @ lat


def wrap_xy(d):
    return d - np.round(d)  # minimal image for in-plane fractional delta


def main():
    if "--selftest" in sys.argv:
        return selftest()
    pos_in, contcar = sys.argv[1], sys.argv[2]
    declared = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else None
    lat0, sp0, cnt0, c0, fl0 = parse(Path(pos_in))
    lat1, sp1, cnt1, c1, fl1 = parse(Path(contcar))
    n = len(sp0)
    if sp0 != sp1 or cnt0 != cnt1:
        raise SystemExit("FAIL: species/counts differ between input and CONTCAR")
    c0c = frac_to_cart(c0, lat0)
    c1c = frac_to_cart(c1, lat1)
    # per-atom Cartesian displacement with minimal in-plane image
    dc = c1c - c0c
    for i in range(n):
        dfrac = wrap_xy((c1[i] - c0[i])[:2])
        dc[i, 0:2] = dfrac @ lat0[0:2, 0:2]
    in_plane = np.abs(dc[:, 0:2]).max(axis=1)
    out_plane = dc[:, 2]

    # 1. symmetry re-audit on the relaxed geometry
    aud = audit(str(Path(contcar)))
    sc = aud.get("symmetry_counts", {})
    n_i, n_c3, n_mv = sc.get("inversion", 0), sc.get("c3", 0), sc.get("mirror_v", 0)

    # 2/3. registry: top-layer rigid in-plane shift vs input
    drift = resid = None
    pt_idx = [i for i, sp in enumerate(sp0) if sp == "Pt"]
    pt_z = sorted({round(float(c0[i][2]), 4) for i in pt_idx})
    if len(pt_z) == 2:
        z0, z1 = pt_z[0], pt_z[1]
        tol = 0.5 * (z1 - z0)

        def layers(cc):
            return ([i for i in range(n) if abs(cc[i, 2] - z0) < tol],
                    [i for i in range(n) if abs(cc[i, 2] - z1) < tol])

        def rigid_shift(cc, btm, top):
            shifts = []
            for j in top:
                cand = [i for i in btm if sp0[i] == sp0[j]]
                best, dmin = None, None
                for i in cand:
                    d = wrap_xy((cc[j][:2] - cc[i][:2]))
                    dist = float(np.linalg.norm(d))
                    if dmin is None or dist < dmin:
                        best, dmin = d, dist
                if best is not None:
                    shifts.append(best)
            return np.median(np.array(shifts), axis=0) if shifts else None

        b0, t0 = layers(c0)
        b1, t1 = layers(c1)
        s0 = rigid_shift(c0, b0, t0)
        s1 = rigid_shift(c1, b1, t1)
        if s0 is not None and s1 is not None:
            drift = float(np.abs(s1 - s0).max())
        if s1 is not None and declared is not None:
            resid = float(np.abs(s1 - np.array(declared)).max())

    # gates: bilayers (x,y FIXED) -> max in-plane displacement < 1e-3 A
    is_bilayer = len(pt_z) == 2
    in_plane_gate = 1e-3
    max_ip = float(in_plane.max()) if is_bilayer else None
    ok = True
    if n_i == 0:
        ok = False
        reason = "inversion lost after relaxation"
    if is_bilayer and max_ip is not None and max_ip > in_plane_gate:
        ok = False
        reason = f"in-plane displacement {max_ip:.2e} A > {in_plane_gate} A"
    if drift is not None and drift > 1e-3:
        ok = False
        reason = f"registry drift {drift:.2e} > 1e-3"
    report = {
        "inversion": bool(n_i), "C3": bool(n_c3), "mirror_v": bool(n_mv),
        "bilayer": is_bilayer,
        "max_in_plane_displacement_A": max_ip,
        "per_atom_in_plane_A": [round(float(x), 6) for x in in_plane],
        "per_atom_out_of_plane_A": [round(float(x), 6) for x in out_plane],
        "registry_drift_frac": drift,
        "slide_after_resid_frac": resid,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"-> {'PASS' if ok else 'FAIL: ' + reason}")
    sys.exit(0 if ok else 1)


def selftest():
    from read_poscar_utils import basic_poscar
    base = Path(__file__).parent / ".." / "INPUTS" / "2253_bilayer_baseline"
    pos = base / "POSCAR"
    # fixture A: INVERSION-PRESERVING z-only relaxation (top +0.02, bottom -0.02
    # in fractional z: mirror-symmetric about the midgap -> inversion retained,
    # in-plane displacements exactly zero)
    lat, sp, cnt, coords, fl = parse(pos)
    c = coords.copy().astype(float)
    z0, z1 = sorted({round(float(c[i][2]), 4) for i in range(len(sp)) if sp[i] == "Pt"})
    tol = 0.5 * (z1 - z0)
    for i in range(len(c)):
        if sp[i] == "Pt":
            c[i, 2] += 0.02 if abs(c[i, 2] - z1) < tol else -0.02
        elif sp[i] == "Se":
            c[i, 2] += 0.01 if abs(c[i, 2] - z1) < tol else -0.01
    fa = Path(__file__).parent / ".." / "LOGS" / "test_contcar_zonly"
    basic_poscar(fa, lat, sp, cnt, c, sd_flags=["F F T"] * len(c))
    main_called = run_verify(str(pos), str(fa), "0", "0")
    assert main_called["PASS"], "fixture A (z-only) should PASS"
    # fixture B: in-plane drift 0.05 A on the top Pt (registry broken)
    c2 = coords.copy().astype(float)
    c2[1, 0] += 0.05  # top Pt in-plane drift (index 1 = top Pt in 2253)
    fb = Path(__file__).parent / ".." / "LOGS" / "test_contcar_drift"
    basic_poscar(fb, lat, sp, cnt, c2)
    run_verify(str(pos), str(fb), "0", "0", expect_fail=True)
    print("verify_relaxed selftest: PASS (z-only PASS; in-plane drift FAIL detected)")
    return 0


def run_verify(pos, cont, ux, uy, expect_fail=False):
    """Re-invoke main logic in-process (avoids exec noise)."""
    old = sys.argv
    sys.argv = ["verify_relaxed", pos, cont, ux, uy]
    res = {"PASS": False}
    try:
        main()
        res["PASS"] = True
    except SystemExit as e:
        res["PASS"] = (e.code == 0)
    sys.argv = old
    if expect_fail and res["PASS"]:
        raise AssertionError("expected FAIL but got PASS")
    if not expect_fail and not res["PASS"]:
        raise AssertionError("expected PASS but got FAIL")
    return res


if __name__ == "__main__":
    main()
