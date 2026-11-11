#!/usr/bin/env python3
"""Stage8k independent nondegenerate spin analysis on raw WAVECAR + EIGENVAL.

Uses stage8k_independent_reader (own parser, not pymatgen). Layout anchored by:
  - file-size closure (97*114528 B per kpoint; verified 4,899,393,312 B)
  - enocc energies vs EIGENVAL per band
  - VASP ncl spinor storage = grouped (up all G, then dn all G) per wavecar_rs/VASP docs;
    interleaved layout reported for contrast only.
Reports pseudo-wavefunction spin (PAW augmentation NOT reconstructible from WAVECAR alone;
PROCAR truncated -> BLOCKED_PAW_SPIN_OPERATOR is a possible hard-stop label).
"""
import json, sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage8k_independent_reader import (read_header, kpoint_meta, band_coeffs,
                                       split_spinor, spin_expected, find_kindex)

NELECT = 36

def parse_eigenval(path):
    with open(path) as f:
        L = f.read().splitlines()
    i = 6
    kpts = []
    bands = []
    for k in range(441):
        while i < len(L) and not L[i].strip():
            i += 1
        kf = [float(v) for v in L[i].split()[:3]]
        i += 1
        rows = []
        for b in range(96):
            v = L[i].split()
            rows.append([float(v[1]), float(v[2])])
            i += 1
        kpts.append(kf)
        bands.append(rows)
    return kpts, bands

TARGETS = {
    "GAMMA": (0, 0, 0),
    "TRUE_K": (1/3, 1/3, 0),
    "TRUE_K_PRIME": (-1/3, -1/3, 0),
    "GAMMA_M_2_3__OLD_FALSE_Kp": (-1/3, 1/3, 0),
    "GAMMA_M_NEAR_M": (0.47619048, 0, 0),
    "GAMMA_K_INTERIOR": (0.14285714, 0.14285714, 0),
    "M": (0.5, 0, 0),
    "LAMBDA": (1/6, 1/6, 0),
}

def analyze(path_w, path_e, label):
    h = read_header(path_w)
    kpts, bands = parse_eigenval(path_e)
    kpts_arr = np.asarray(kpts)
    res = {"job": label, "header": {k: h[k] for k in ("recl", "spin", "rtag", "nk", "nb", "encut", "efermi")},
           "parser": "stage8k_independent_reader (raw f64/complex64 record walk; grouped ncl spinor)",
           "per_target": {}}
    for name, tgt in TARGETS.items():
        ki, dist = find_kindex(h, tgt, kpts_arr)
        nplane, wk, enocc = kpoint_meta(path_w, h, ki)
        E_wc = enocc[:, 0]
        E_eig = np.array([bands[ki][b][0] for b in range(h["nb"])])
        dE = E_wc - E_eig
        entry = {"k_index": ki, "WAVECAR_k": [float(x) for x in wk],
                 "EIGENVAL_k": kpts[ki], "dist": round(float(dist), 6),
                 "nplane": nplane,
                 "is_exact_requested_k_on_mesh": bool(dist < 1e-5),
                 "k_label_warning": ("nearest grid point, NOT the exact requested high-symmetry coordinate" if dist >= 1e-5 else None),
                 "max_enocc_energy_dev_eV": round(float(np.max(np.abs(dE))), 6),
                 "energy_crosscheck": {},
                 "bands_near_gap": []}
        for b in (33, 34, 35, 36, 37):
            entry["energy_crosscheck"][f"band{b+1}"] = [
                round(float(E_wc[b]), 6), round(float(E_eig[b]), 6)]
        # band isolation from FULL EIGENVAL spectrum at this k
        for b in (33, 34, 35, 36, 37):
            e = E_eig[b]
            below = E_eig[b-1] if b > 0 else None
            above = E_eig[b+1] if b + 1 < h["nb"] else None
            gap = min([x for x in (None if below is None else abs(e-below),
                                   None if above is None else abs(e-above)) if x is not None])
            entry["bands_near_gap"].append({
                "band": b + 1, "E_eV": round(float(e), 6),
                "min_adjacent_EIGENVAL_gap_eV": round(float(gap), 6),
                "occ": round(float(bands[ki][b][1]), 4)})
        # spinor expectations: grouped (physical) and interleaved (contrast)
        spin_rows = []
        for b in (33, 34, 35, 36, 37):
            coeffs, _ = band_coeffs(path_w, h, ki, b)
            row = {"band": b + 1, "E_eV": round(float(E_eig[b]), 6)}
            for layout in ("grouped", "interleaved"):
                up, dn, nG = split_spinor(coeffs, nplane, layout)
                sx, sy, sz, tot = spin_expected(up, dn)
                row[f"S_{layout}"] = ([round(sx, 6), round(sy, 6), round(sz, 6)]
                                      if sx is not None else None)
                row[f"norm_{layout}"] = round(float(tot), 6)
            # phase invariance on band35
            if b == 34:
                up, dn, _ = split_spinor(coeffs, nplane, "grouped")
                sx0, sy0, sz0, _ = spin_expected(up, dn)
                phases = []
                for phi in (0.7, 1.9, 3.1):
                    up2 = up * np.exp(1j * phi)
                    dn2 = dn * np.exp(1j * phi)
                    sx, sy, sz, _ = spin_expected(up2, dn2)
                    phases.append([round(float(sx), 6), round(float(sy), 6), round(float(sz), 6)])
                row["phase_invariance_band35"] = {
                    "base": [round(sx0, 6), round(sy0, 6), round(sz0, 6)], "phases": phases}
            spin_rows.append(row)
        entry["spinor_rows"] = spin_rows
        res["per_target"][name] = entry
    # time reversal: TRUE_K vs TRUE_K_PRIME spin of band35/36 (grouped)
    tr = {}
    for b in (35, 36):
        kiK, _ = find_kindex(h, TARGETS["TRUE_K"], kpts_arr)
        kiKp, _ = find_kindex(h, TARGETS["TRUE_K_PRIME"], kpts_arr)
        cK, npK = band_coeffs(path_w, h, kiK, b - 1)
        cKp, npKp = band_coeffs(path_w, h, kiKp, b - 1)
        upK, dnK, _ = split_spinor(cK, npK, "grouped")
        upKp, dnKp, _ = split_spinor(cKp, npKp, "grouped")
        sK = spin_expected(upK, dnK)
        sKp = spin_expected(upKp, dnKp)
        tr[f"band{b}"] = {"K_S": [round(x, 6) if x is not None else None for x in sK[:3]],
                          "Kp_S": [round(x, 6) if x is not None else None for x in sKp[:3]],
                          "E_K": round(float(E_eig[0] * 0 + bands[kiK][b-1][0]), 6),
                          "E_Kp": round(float(bands[kiKp][b-1][0]), 6),
                          "note": "nonmagnetic TR: E(K)=E(-K) required; spin vectors need not equal -S without band pairing proof"}
    res["time_reversal_K_vs_Kp"] = tr
    return res

if __name__ == "__main__":
    import os
    root = os.environ.get("STAGE8I_ROOT", "/work/sunq/xjq/runs-stage8i-20261021")
    out = {}
    for lab, sub in (("soc_AB_s0.00", "soc_AB_s0.00"), ("soc_GLIDE_s0.50", "soc_GLIDE_s0.50")):
        print("analyzing", sub, flush=True)
        out[lab] = analyze(f"{root}/{sub}/WAVECAR", f"{root}/{sub}/EIGENVAL", lab)
    out_path = os.environ.get("STAGE8K_OUT", "/work/sunq/xjq/stage8i-work/tools/stage8k_independent_spin.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    print("WROTE", out_path)
