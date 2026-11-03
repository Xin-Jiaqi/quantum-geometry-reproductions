#!/usr/bin/env python3
"""parse_vasp_b1_cg.py — read-only forensic parser for the AgF2 B1 WSL CG relaxation.

Reads VASP 6.6.0 raw outputs (OUTCAR / OSZICAR / vasprun.xml / stdout) and the
starting POSCAR, and emits:
  B1_CG_TRAJECTORY.csv   — per-ionic-step quantitative table
  B1_CG_ION_FORCES.csv   — 15x12 (atom, step) forces & coordinates
  B1_CG_AUDIT.json       — cross-checks and evidence-level tagging

Only Python stdlib is used. No VASP binary, no POTCAR content is read.
Atom ordering (1-based): Ag 1-4 (top 1,2 / bottom 3,4); F 5-12
(top-half F 5-8: inner 5,8 outer 6,7; bottom-half F 9-12: inner 9,12 outer 10,11).

Usage:
  python3 parse_vasp_b1_cg.py OUTCAR OSZICAR vasprun.xml vasp.stdout POSCAR [--outdir DIR]

Test:
  python3 parse_vasp_b1_cg.py --selftest
"""
import argparse
import json
import math
import os
import re
import sys
from collections import OrderedDict

# ----------------------------------------------------------------------------
# VASP file parsing (version-tolerant: missing fields -> None, never invented)
# ----------------------------------------------------------------------------

class IonStep:
    """One completed ionic step (1-based index), with raw values as found."""
    def __init__(self, idx):
        self.idx = idx
        self.positions = None      # 12 x 3 cartesian (A)
        self.forces = None         # 12 x 3 (eV/A)
        self.toten = None
        self.energy_no_entropy = None
        self.energy_sigma0 = None
        self.ediff_reached = None
        self.oszicar = None        # dict F/E0/dE/mag
        self.total_magnetization = None
        self.vasprun_forces = None # 12 x 3 (eV/A) from vasprun.xml
        self.vasprun_e0 = None

def parse_outcar(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    steps = []
    # force blocks
    blocks = re.findall(r"TOTAL-FORCE \(eV/Angst\)\n.*?\n(.*?)\n\s*total drift", txt, re.S)
    for b in blocks:
        rows = []
        for l in b.splitlines():
            p = l.split()
            if len(p) == 6:
                try:
                    v = [float(x) for x in p]
                    rows.append((v[:3], v[3:]))
                except ValueError:
                    pass
        steps.append(rows)
    # energies per step
    toten = [float(x) for x in re.findall(r"free  energy   TOTEN\s*=\s*([-\d.]+)", txt)]
    eno   = re.findall(r"energy  without entropy=\s+([-\d.]+)\s+energy\(sigma->0\)\s*=\s+([-\d.]+)", txt)
    # EDIFF markers are not per-step anchored; count them for the whole run
    n_ediff = txt.count("aborting loop because EDIFF is reached")
    # magnetization (x): take the LAST occurrence block only (final step's
    # per-atom moments); per-step net M comes from OSZICAR. rfind+slice avoids
    # catastrophic backtracking on large files.
    mags = []
    mi = txt.rfind("magnetization (x)")
    if mi != -1:
        seg = txt[mi:mi + 4000]
        lines = seg.splitlines()
        tot = {}
        netm = None
        for l in lines[2:]:
            p = l.split()
            if len(p) >= 5 and p[0].replace(".", "").replace("-", "").isdigit():
                try:
                    tot[int(p[0])] = float(p[1])
                except ValueError:
                    pass
            elif l.startswith("   total magnetization"):
                try:
                    netm = float(l.split()[-1])
                except ValueError:
                    pass
        mags.append((tot, netm))
    nsteps = len(steps)
    for i in range(nsteps):
        s = IonStep(i + 1)
        s.positions = [r[0] for r in steps[i]]
        s.forces = [r[1] for r in steps[i]]
        s.toten = toten[i] if i < len(toten) else None
        if i < len(eno):
            s.energy_no_entropy = float(eno[i][0])
            s.energy_sigma0 = float(eno[i][1])
        if i < len(mags):
            s.total_magnetization = mags[0][1]
        s.ediff_reached = n_ediff >= nsteps  # whole-run count >= step count
        steps[i] = s
    return steps, n_ediff

def parse_oszicar(path):
    out = OrderedDict()
    pat = re.compile(
        r"^\s*(\d+)\s+F=\s*([-\d.Ee+]+)\s+E0=\s*([-\d.Ee+]+)\s+d E\s*=\s*([-\d.Ee+]+)\s+mag=\s*([-\d.Ee+]+)")
    for l in open(path, encoding="utf-8", errors="replace"):
        m = pat.match(l)
        if m:
            out[int(m.group(1))] = {
                "F": float(m.group(2)), "E0": float(m.group(3)),
                "dE": float(m.group(4)), "mag": float(m.group(5)),
            }
    return out

def parse_vasprun(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    # forces per ionic step: <varray name="forces"> blocks; linear scan via
    # finditer + bounded slices avoids catastrophic backtracking on big files.
    forces = []
    for m in re.finditer(r"<varray name=\"forces\"[^>]*>", txt):
        start = m.end()
        end = txt.find("</varray>", start)
        rows = []
        for v in re.finditer(r"<v>\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*</v>", txt[start:end]):
            rows.append([float(v.group(1)), float(v.group(2)), float(v.group(3))])
        forces.append(rows)
    # VASP 6.x: e_fr_energy is written per <scstep>; the LAST scstep of each
    # <calculation> equals that ionic step's TOTEN (free energy).
    segs = txt.split("<calculation>")
    fr = []
    for seg in segs[1:]:
        vals = re.findall(r"<i name=\"e_fr_energy\">\s*([-\d.Ee+]+)\s*</i>", seg)
        if vals:
            fr.append(float(vals[-1]))
    return {"forces": forces, "e0": fr}

def parse_stdout_markers(path):
    """Extract CG line-search markers with stdout line numbers."""
    markers = []
    pats = [
        r"^\s*(curvature:[^\n]*)",
        r"^\s*(trial:[^\n]*)",
        r"^\s*(search vector[^\n]*)",
        r"^\s*(ZBRENT:[^\n]*)",
        r"^\s*(opt\s*:[^\n]*)",
        r"^\s*(step\s*:[^\n]*)",
        r"^\s*(trial-energy change:[^\n]*)",
        r"^\s*(soft stop[^\n]*)",
        r"^\s*(deleting file STOPCAR[^\n]*)",
    ]
    for ln, l in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
        for pat in pats:
            m = re.match(pat, l)
            if m:
                markers.append({"line": ln, "text": m.group(1).strip()})
                break
    return markers

# ----------------------------------------------------------------------------
# Geometry helpers (orthorhombic Ag4F8 bilayer; PBC min-image on x,y; no z wrap)
# ----------------------------------------------------------------------------

def geometry_metrics(pos, lat, ref=None):
    a, b, c = lat
    def layer_mean(ids):
        return [sum(pos[i - 1][k] for i in ids) / len(ids) for k in range(3)]
    def minimg_disp(p, q):
        d = [p[k] - q[k] for k in range(3)]
        d[0] -= a * round(d[0] / a)
        d[1] -= b * round(d[1] / b)
        d[2] -= c * round(d[2] / c)
        return d
    def norm3(d):
        return math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    top_ag = layer_mean([1, 2])
    bot_ag = layer_mean([3, 4])
    top_inner_f = layer_mean([5, 8])
    bot_inner_f = layer_mean([9, 12])
    ff = min(norm3(minimg_disp(pos[i - 1], pos[j - 1]))
             for i in (5, 6, 7, 8) for j in (9, 10, 11, 12))
    disp_y = top_ag[1] - bot_ag[1]
    disp_y -= b * round(disp_y / b)
    # interlayer relative sliding: per-atom y displacement vs the reference
    # (step-1) geometry, unwrapped by +/- b to stay nearest, then layer mean
    # difference. This avoids mixing the two Ag sublattice images.
    slide = None
    if ref is not None:
        def dy_mean(ids):
            ds = []
            for i in ids:
                d = pos[i - 1][1] - ref[i - 1][1]
                d -= b * round(d / b)
                ds.append(d)
            return sum(ds) / len(ds)
        slide = dy_mean([1, 2]) - dy_mean([3, 4])
    return {
        "ag_z_sep_A": top_ag[2] - bot_ag[2],
        "top_ag_z_A": top_ag[2],
        "bot_ag_z_A": bot_ag[2],
        "inner_f_gap_z_A": top_inner_f[2] - bot_inner_f[2],
        "nearest_interlayer_f_f_A": ff,
        "ag_y_registry_A": disp_y,
        "ag_y_registry_frac_b": disp_y / b,
        "ag_y_relative_slide_A": slide,
    }

def parse_poscar(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    lat = [[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    mode = lines[7].strip().lower()
    if mode.startswith("d"):
        frac = [[float(x) for x in l.split()[:3]] for l in lines[8:8 + n]]
        cart = [[frac[i][k] * lat[k][k] for k in range(3)] for i in range(n)]
    else:
        cart = [[float(x) for x in l.split()[:3]] for l in lines[8:8 + n]]
    return {"lat": [lat[0][0], lat[1][1], lat[2][2]], "cart": cart, "counts": counts}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def run(outcar, oszicar, vasprun, stdout, poscar, outdir=None):
    steps, n_ediff = parse_outcar(outcar)
    osz = parse_oszicar(oszicar)
    vr = parse_vasprun(vasprun)
    markers = parse_stdout_markers(stdout)
    start = parse_poscar(poscar)
    a, b, c = start["lat"]

    # align vasprun forces: first <varray name="forces"> appears for the final
    # structure of the first SCF (step 1) — take the LAST nsteps entries
    vr_forces = vr["forces"][-len(steps):] if len(vr["forces"]) >= len(steps) else None
    max_abs_force_diff = None
    rows = []
    for s in steps:
        m = geometry_metrics(s.positions, (a, b, c), ref=start["cart"])
        fmax = max(s.forces, key=lambda f: math.sqrt(f[0] ** 2 + f[1] ** 2 + f[2] ** 2))
        fmax_norm = math.sqrt(sum(x * x for x in fmax))
        amax = s.forces.index(fmax) + 1
        row = OrderedDict([
            ("step", s.idx),
            ("toten_ev", s.toten),
            ("e0_ev", s.energy_no_entropy),
            ("e_sigma0_ev", s.energy_sigma0),
            ("oszicar_F_ev", osz.get(s.idx, {}).get("F")),
            ("oszicar_E0_ev", osz.get(s.idx, {}).get("E0")),
            ("oszicar_dE_ev", osz.get(s.idx, {}).get("dE")),
            ("oszicar_mag_mub", osz.get(s.idx, {}).get("mag")),
            ("electronic_ediFF_converged", bool(s.ediff_reached)),
            ("max_force_norm_ev_ang", round(fmax_norm, 4)),
            ("max_abs_Fz_ev_ang", round(max(abs(f[2]) for f in s.forces), 4)),
            ("max_abs_Fy_ev_ang", round(max(abs(f[1]) for f in s.forces), 4)),
            ("max_force_atom_1based", amax),
            ("ag_z_sep_A", round(m["ag_z_sep_A"], 4)),
            ("inner_f_gap_z_A", round(m["inner_f_gap_z_A"], 4)),
            ("nearest_interlayer_f_f_A", round(m["nearest_interlayer_f_f_A"], 4)),
            ("ag_y_registry_A", round(m["ag_y_registry_A"], 4)),
            ("ag_y_registry_frac_b", round(m["ag_y_registry_frac_b"], 6)),
            ("ag_y_relative_slide_A", round(m["ag_y_relative_slide_A"], 6) if m["ag_y_relative_slide_A"] is not None else None),
            ("total_magnetization_mub", s.total_magnetization),
            ("vasprun_e0_ev", vr["e0"][-len(steps) + s.idx - 1] if len(vr["e0"]) >= len(steps) else None),
        ])
        if vr_forces and s.idx <= len(vr_forces):
            diffs = [math.sqrt(sum((vr_forces[s.idx - 1][i][k] - s.forces[i][k]) ** 2 for k in range(3)))
                     for i in range(12)]
            row["outcar_vs_vasprun_max_abs_force_diff_ev_ang"] = round(max(diffs), 6)
            if max_abs_force_diff is None or max(diffs) > max_abs_force_diff:
                max_abs_force_diff = max(diffs)
        rows.append(row)

    # ion forces CSV
    ion_rows = []
    for s in steps:
        for i in range(12):
            ion_rows.append({
                "step": s.idx, "atom_1based": i + 1,
                "species": "Ag" if i < 4 else "F",
                "x_A": s.positions[i][0], "y_A": s.positions[i][1], "z_A": s.positions[i][2],
                "Fx_ev_ang": s.forces[i][0], "Fy_ev_ang": s.forces[i][1], "Fz_ev_ang": s.forces[i][2],
            })

    if outdir:
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, "B1_CG_TRAJECTORY.csv"), "w") as f:
            f.write(",".join(rows[0].keys()) + "\n")
            for r in rows:
                f.write(",".join("" if v is None else str(v) for v in r.values()) + "\n")
        with open(os.path.join(outdir, "B1_CG_ION_FORCES.csv"), "w") as f:
            f.write(",".join(ion_rows[0].keys()) + "\n")
            for r in ion_rows:
                f.write(",".join(str(r[k]) for k in ion_rows[0].keys()) + "\n")

    return {
        "n_ionic_steps": len(steps),
        "n_ediFF_markers_total": n_ediff,
        "n_oszicar_steps": len(osz),
        "n_vasprun_force_blocks": len(vr["forces"]),
        "n_stdout_markers": len(markers),
        "trajectory": rows,
        "ion_forces": ion_rows,
        "markers": markers,
        "outcar_vs_vasprun_max_abs_force_diff_ev_ang": max_abs_force_diff,
    }

def selftest():
    """Synthetic minimal fixtures validate parsing + geometry without VASP."""
    import tempfile
    d = tempfile.mkdtemp()
    # POSCAR-like: 2-layer toy cell with 12 atoms (Ag4F8 ordering, orthorhombic)
    lat = (4.0, 5.0, 20.0)
    cart = [[0.0, 0.0, 15.0], [2.0, 0.0, 15.0], [0.0, 2.5, 10.0], [2.0, 2.5, 10.0],
            [1.0, 1.0, 14.0], [3.0, 1.0, 16.0], [1.0, 4.0, 16.0], [3.0, 4.0, 14.0],
            [1.0, 1.0, 11.0], [3.0, 1.0, 9.0], [1.0, 4.0, 9.0], [3.0, 4.0, 11.0]]
    pos = ["toy", "1.0", "%.6f 0 0" % lat[0], "0 %.6f 0" % lat[1], "0 0 %.6f" % lat[2],
           "Ag F", "4 8", "Cartesian"] + ["%.6f %.6f %.6f" % tuple(p) for p in cart]
    ppos = os.path.join(d, "POSCAR")
    open(ppos, "w").write("\n".join(pos) + "\n")
    # OUTCAR: 2 force blocks + energies
    out = ["toy outcar"]
    for i in range(2):
        out.append("TOTAL-FORCE (eV/Angst)")
        out.append("-" * 70)
        for j in range(12):
            out.append("%.6f %.6f %.6f %.6f %.6f %.6f" % (cart[j][0], cart[j][1], cart[j][2],
                                                        0.1 * (j + 1), 0.0, 0.01 * (j + 1)))
        out.append("   total drift: 0.0 0.0 0.0")
        out.append("free  energy   TOTEN  = -%d.1234 eV" % (10 + i))
        out.append("energy  without entropy= -%d.1234  energy(sigma->0) = -%d.1234" % (10 + i, 10 + i))
        out.append("magnetization (x)")
        out.append("   total magnetization = %.4f" % (0.5 * (i + 1)))
    out.append("aborting loop because EDIFF is reached")
    pout = os.path.join(d, "OUTCAR")
    open(pout, "w").write("\n".join(out) + "\n")
    # OSZICAR
    poz = os.path.join(d, "OSZICAR")
    open(poz, "w").write("\n".join(
        "   %d F= -%d.1234 E0= -%d.1234  d E =-%d.1234  mag=  0.0000" % (i + 1, 10 + i, 10 + i, i)
        for i in range(2)) + "\n")
    # vasprun.xml: 2 force varrays + e_0
    pv = os.path.join(d, "vasprun.xml")
    with open(pv, "w") as f:
        f.write("<vasprun><i name='e_0'>-11.0</i><i name='e_0'>-12.0</i>")
        for i in range(2):
            f.write("<varray name='forces'>")
            for j in range(12):
                f.write("<v>%f %f %f</v>" % (0.1 * (j + 1), 0.0, 0.01 * (j + 1)))
            f.write("</varray>")
        f.write("</vasprun>")
    # stdout markers
    ps = os.path.join(d, "vasp.stdout")
    open(ps, "w").write("  1 F= -10.1234\ncurvature:  -1.00\nZBRENT: bracketing found\n"
                        "opt :  3.0000  next Energy=   -11.000\nsoft stop encountered!\n")
    res = run(pout, poz, pv, ps, ppos, outdir=d)
    assert res["n_ionic_steps"] == 2, res
    assert res["trajectory"][0]["max_force_norm_ev_ang"] == round(math.sqrt(1.2 ** 2 + 0.12 ** 2), 4)
    assert res["trajectory"][0]["max_force_atom_1based"] == 12
    assert abs(res["trajectory"][0]["ag_z_sep_A"] - (15.0 - 10.0)) < 1e-9
    assert abs(res["trajectory"][0]["inner_f_gap_z_A"] - (14.0 - 11.0)) < 1e-9
    print("SELFTEST OK — parsed %d steps, %d markers" % (res["n_ionic_steps"], len(res["markers"])))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="OUTCAR OSZICAR vasprun.xml vasp.stdout POSCAR")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if len(args.files) != 5:
        ap.error("need exactly 5 files: OUTCAR OSZICAR vasprun.xml vasp.stdout POSCAR")
    res = run(*args.files, outdir=args.outdir)
    print(json.dumps({
        "n_ionic_steps": res["n_ionic_steps"],
        "n_ediFF_markers_total": res["n_ediFF_markers_total"],
        "n_oszicar_steps": res["n_oszicar_steps"],
        "n_vasprun_force_blocks": res["n_vasprun_force_blocks"],
        "n_stdout_markers": len(res["markers"]),
        "outcar_vs_vasprun_max_abs_force_diff_ev_ang": res["outcar_vs_vasprun_max_abs_force_diff_ev_ang"],
    }, indent=1))

if __name__ == "__main__":
    main()
