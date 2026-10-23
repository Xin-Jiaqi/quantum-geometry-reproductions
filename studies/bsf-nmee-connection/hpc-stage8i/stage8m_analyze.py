#!/usr/bin/env python3
"""Stage8m post-analysis: convergence and symmetry checks over the response
matrix (window x B-step x lambda-step x source x axis). Pure JSON processing.

Checks:
 1. lambda-step stability: for fixed (src, win, B), Gamma_BZ across dl values
 2. B-step (central-difference) convergence: |Gamma(B2) - Gamma(B1)| trend
 3. nested-window / remote-band truncation stability across 36..96
 4. glide total-response null: |Gamma| of glide vs AB at matching configs
 5. in-plane tensor symmetry residual |Gamma_x + Gamma_y|
Usage: stage8m_analyze.py --ab ab.json --gl glide.json [--out report.json]
"""
import argparse, json, math


def gval(cfg, axis):
    return cfg["axis_" + axis]["Gamma_BZ_Ang2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ab", required=True)
    ap.add_argument("--gl", required=True)
    ap.add_argument("--out", default="stage8m_analysis.json")
    a = ap.parse_args()
    ab = json.load(open(a.ab))
    gl = json.load(open(a.gl))
    ab_cfg = {k: v for k, v in ab.items() if k != "_meta"}
    gl_cfg = {k: v for k, v in gl.items() if k != "_meta"}
    # parse config keys: src{x,y}_win{W}_B{B}_dl{D}
    def parse(k):
        parts = {}
        for tok in k.split("_"):
            if tok.startswith("src"):
                parts["src"] = tok[3:]
            elif tok.startswith("win"):
                parts["win"] = int(tok[3:])
            elif tok.startswith("B"):
                parts["B"] = float(tok[1:])
            elif tok.startswith("dl"):
                parts["dl"] = float(tok[2:])
        return parts

    keys = sorted(k for k in ab_cfg if not k.startswith("_"))
    rows = []
    # B convergence: for each (src, win, dl): ratio G(B=.25)/G(B=.125)
    bconv = {}
    dlconv = {}
    for k in keys:
        p = parse(k)
        if p["B"] == 0.125:
            k2 = k.replace("B0.125", "B0.25")
            if k2 in ab_cfg:
                for ax in ("x", "y"):
                    g1, g2 = gval(ab_cfg[k], ax), gval(ab_cfg[k2], ax)
                    bconv[(p["src"], p["win"], p["dl"], ax)] = (g1, g2)
    # dl convergence: for fixed (src, win, B): values across dl
    for k in keys:
        p = parse(k)
        if p["dl"] == 1e-4:
            dlrows = []
            for dl in (1e-4, 2e-4, 4e-4):
                kk = k.replace("dl0.0001", f"dl{dl:g}")
                if kk in ab_cfg:
                    dlrows.append((dl, gval(ab_cfg[kk], "x"),
                                   gval(ab_cfg[kk], "y")))
            if len(dlrows) == 3:
                dlconv[(p["src"], p["win"], p["B"])] = dlrows
    # glide null: matched configs
    gl_null = {}
    for k in keys:
        if k in gl_cfg:
            for ax in ("x", "y"):
                gl_null[(k, ax)] = (gval(ab_cfg[k], ax), gval(gl_cfg[k], ax))
    # window scan at reference config
    ref = {}
    for src in ("x", "y"):
        ref[src] = []
        for win in (36, 48, 64, 80, 96):
            k = f"src{src}_win{win}_B0.25_dl0.0001"
            if k in ab_cfg:
                ref[src].append((win, gval(ab_cfg[k], "x"),
                                 gval(ab_cfg[k], "y")))
    # symmetry residual per config
    sym = {}
    for k in keys:
        sym[k] = ab_cfg[k]["sym_residual_abs"]
    max_sym = max(sym.values()) if sym else float("nan")
    report = {
        "configs": len(keys),
        "b_step_convergence": {
            "rows": [{"key": k, "Gx_B125": g1, "Gx_B25": g2}
                     for (src, win, dl, ax), (g1, g2) in
                     sorted(bconv.items()) if ax == "x"][:40],
            "note": ("B-step is a central difference at fixed lambda; the "
                     "ratio G(.25)/G(.125) should be ~1 if the derivative is "
                     "converged (it is a field derivative, not a power)."),
        },
        "lambda_step_stability": {
            "rows": [{"key": k, "dl": d, "Gx": gx, "Gy": gy}
                     for k, dlrows in sorted(dlconv.items())
                     for d, gx, gy in dlrows][:60],
            "note": "Gamma/dl should be flat across the dl sequence.",
        },
        "window_scan_at_B0.25_dl1e-4": {
            "src_x": ref.get("x", []),
            "src_y": ref.get("y", []),
            "note": ("Nested windows 36/48/64/80/96 test remote-band "
                     "truncation; a stable tail indicates the response is "
                     "not carried by distant bands."),
        },
        "glide_total_null": {
            "max_abs_glide_Gx": max((abs(v) for k, ax in gl_null
                                     for v in (gl_null[(k, ax)][1],) if ax == "x"),
                                    default=float("nan")),
            "max_abs_AB_Gx": max((abs(gl_null[(k0, "x")][0])
                                  for (k0, ax) in gl_null if ax == "x"),
                                 default=float("nan")),
            "rows": [{"key": k, "AB": gl_null[(k, ax)][0],
                      "GL": gl_null[(k, ax)][1]}
                     for k, ax in sorted(gl_null) if ax == "x"][:30],
            "note": "layer-exchanging glide should give total response 0.",
        },
        "sym_residual_max_abs": max_sym,
    }
    json.dump(report, open(a.out, "w"), indent=1)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
