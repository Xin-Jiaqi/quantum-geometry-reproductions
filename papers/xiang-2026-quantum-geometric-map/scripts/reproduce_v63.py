"""Reproduce v63 non-invasive QWZ voltage-probe and finite-size controls."""
import argparse
import json
from pathlib import Path
import numpy as np
from qgmap.qwz_probe_limit_v63 import probe_strength_response
from qgmap.qwz_charge_kubo_v62 import qwz_charge_kubo


def compact(d, include_cuts=False):
    out = {
        "L": d["L"], "drive_width": d["width"], "probe_width": d["probe_width"],
        "eta": d["eta"], "probe_scale": d["probe_scale"],
        "Rxy_h_over_q2": d["Rxy_h_over_q2"],
        "drive_current_q2_over_h": d["drive_current_q2_over_h"],
        "top_voltage": float(d["voltages"][1]), "bottom_voltage": float(d["voltages"][3]),
        "max_cut_deviation": float(np.max(np.abs(d["x_cut_currents"]-d["drive_current_q2_over_h"]))),
        "top_probe_gamma_norm": float(d["gamma_norms"][1]),
    }
    if include_cuts:
        out["x_cut_currents"] = [float(x) for x in d["x_cut_currents"]]
    return out


def create():
    out = {
        "version": "v63",
        "scientific_status": "VERIFIED specified finite-device identities; CONDITIONAL non-invasive probe limit; UNRESOLVED generic Kubo-NEGF limit and original orbital mixed E/B",
        "model": "same QWZ finite device and semi-infinite QWZ leads; top/bottom coupling amplitude lambda; Sigma/Gamma scale lambda^2",
        "temperature": "zero-temperature coherent finite-device Fermi-energy response",
        "convention": "V_L-V_R=1; I_T=I_B=0; Rxy=(V_top-V_bottom)/I_left",
        "units": {"Rxy": "h/q^2", "conductance": "q^2/h per unit voltage", "eta": "dimensionless QWZ energy"},
        "strength_sweep": {},
        "profile_examples": {},
        "probe_width_sweep_at_unit_scale": {},
        "probe_width_sweep_weak_scale": {},
        "size_sweep_weak_scale": {},
        "eta_sweep_weak_scale": {},
        "position_control": {},
        "metal_negative_control": {},
    }
    for scale in (1., .3, .1, .03, .01, .003, .001, .0001):
        d = probe_strength_response(7, -1, 3, eta=1e-4, probe_scale=scale, probe_width=3)
        row = compact(d)
        out["strength_sweep"][str(scale)] = {k: row[k] for k in (
            "probe_scale", "Rxy_h_over_q2", "drive_current_q2_over_h",
            "max_cut_deviation", "top_probe_gamma_norm")}
        if scale in (1., .0001):
            out["profile_examples"][str(scale)] = [float(x) for x in d["x_cut_currents"]]
    for scale, dest in ((1., "probe_width_sweep_at_unit_scale"), (1e-4, "probe_width_sweep_weak_scale")):
        for width in (2, 3, 5, 7):
            d = probe_strength_response(7, -1, 3, eta=1e-4, probe_scale=scale, probe_width=width)
            row = compact(d)
            out[dest][str(width)] = {k: row[k] for k in (
                "Rxy_h_over_q2", "drive_current_q2_over_h", "max_cut_deviation")}
    for L in (5, 7, 9, 11):
        d = probe_strength_response(L, -1, 3, eta=1e-4, probe_scale=1e-4, probe_width=3)
        row = compact(d)
        out["size_sweep_weak_scale"][str(L)] = {k: row[k] for k in (
            "Rxy_h_over_q2", "drive_current_q2_over_h")}
    for eta in (1e-3, 1e-4, 1e-5, 1e-6):
        d = probe_strength_response(7, -1, 3, eta=eta, probe_scale=1e-3, probe_width=3)
        out["eta_sweep_weak_scale"][str(eta)] = d["Rxy_h_over_q2"]
    for offsets in ((0, 0, 0, 0), (0, 1, 0, -1)):
        d = probe_strength_response(7, -1, 3, eta=1e-4, probe_scale=1e-4,
                                    probe_width=3, offsets=offsets)
        out["position_control"][str(offsets)] = d["Rxy_h_over_q2"]
    d = probe_strength_response(5, -1, 3, mu=1.5, eta=1e-4,
                                probe_scale=1e-4, probe_width=3)
    sigma = qwz_charge_kubo(-1, 1.5, 0., 31)["sigma_interband_q2_over_h"]
    out["metal_negative_control"] = {
        "device_Rxy_h_over_q2": d["Rxy_h_over_q2"],
        "bulk_sigma_xy_q2_over_h": sigma,
        "reciprocal_bulk_sigma_h_over_q2": 1/sigma,
        "explicitly_not_claimed_equal": True,
    }
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", default="results/frozen/v63_noninvasive_probe_limit.json")
    args = p.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(create(), sort_keys=True, indent=2)+"\n")
    print(f"v63 frozen evidence -> {path}")
