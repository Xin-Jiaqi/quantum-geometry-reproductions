"""Independent stage-3 gate: electronic band-sum landscape over a full sliding period.

NOT a DFT total-energy, phonon, elastic or ferroelectric stability calculation.
A one-period scan can *disqualify* a claim that given u are observed
minima of this electronic-only toy model; it cannot validate real BSF.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from bsf_igmc_micro.model import (
    Parameters, band_energy, geometric_response, layer_polarization,
)


def landscape(*, mesh: int = 13, shell: int = 3, n_intervals: int = 40):
    """Return data + audit diagnostics. Never treat band sums as total energies."""
    p = replace(Parameters(mass=0.85), shell=shell)
    u_grid = np.linspace(0.0, 1.0, n_intervals+1)
    records = []
    checks = {}
    for branch, mu in (("AF", 0.0), ("FM", -0.03)):
        energy = np.array([band_energy(float(u), branch, p, N=mesh) for u in u_grid])
        for u, e in zip(u_grid, energy):
            records.append({"branch": branch, "u": float(u), "electronic_band_sum_eV_cell": float(e)})
        # Report only properties of this sampled electronic-only model.
        center = len(u_grid)//2
        positive_side = np.diff(energy[:center+1])
        negative_side = np.diff(energy[center:])
        at_u = {}
        for u in (0.0, .19, .5, .81, 1.0):
            ge = geometric_response(u, branch, p, N=mesh, mu=mu)
            at_u[str(u)] = {
                "Dxx_over_omega": ge["Dxx_over_omega"],
                "indirect_gap_eV": ge["indirect_gap"],
                "layer_tau_z": layer_polarization(u, branch, p, N=mesh, mu=mu),
            }
        checks[branch] = {
            "grid_min_u": float(u_grid[np.argmin(energy)]),
            "grid_min_electronic_band_sum_eV_cell": float(np.min(energy)),
            "E_at_u0_eV_cell": float(energy[0]),
            "E_at_u0p5_eV_cell": float(energy[center]),
            "E_at_u1_eV_cell": float(energy[-1]),
            "Delta_E_0p19_minus_0p5_eV_cell": float(
                band_energy(.19, branch, p, N=mesh)-energy[center]
            ),
            "max_periodicity_residual_eV_cell": float(abs(energy[0]-energy[-1])),
            "all_steps_0_to_0p5_decrease": bool(np.all(positive_side < 0)),
            "all_steps_0p5_to_1_increase": bool(np.all(negative_side > 0)),
            "points": at_u,
        }
    return {
        "status": "ELECTRONIC_BAND_SUM_ONLY_NO_FE_CLAIM",
        "parameter_scope": "8-band two-center hopping, no onsite layer polarity, N mesh midpoint",
        "mesh": mesh,
        "shell": shell,
        "n_intervals": n_intervals,
        "grid_0_to_1": True,
        "note": "No ionic/electrostatic/elastic/relaxation energy; u=0.5 and -0.5 are translation equivalent; only a sampled 1D translation path, not 2D landscape.",
        "checks": checks,
    }, records


def main() -> None:
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--mesh", type=int, default=13)
    cli.add_argument("--shell", type=int, default=3)
    cli.add_argument("--intervals", type=int, default=40)
    cli.add_argument("--outdir", type=Path, default=Path("results"))
    args = cli.parse_args()
    if args.mesh < 9 or args.mesh % 2 == 0 or args.shell < 2 or args.intervals < 8 or args.intervals % 2:
        cli.error("Require odd mesh>=9, shell>=2 and even intervals>=8")
    report,records=landscape(mesh=args.mesh,shell=args.shell,n_intervals=args.intervals)
    args.outdir.mkdir(parents=True,exist_ok=True)
    json_path=args.outdir/'stage3_electronic_energy_gate.json'
    csv_path=args.outdir/'stage3_energy_scan_one_period.csv'
    json_path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    with csv_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['branch','u','electronic_band_sum_eV_cell'])
        writer.writeheader();writer.writerows(records)
    print(json.dumps({'result':str(json_path),'csv':str(csv_path),'checks':{
        b:{k:v for k,v in a.items() if k!='points'} for b,a in report['checks'].items()},
        'status':report['status']},indent=2))

if __name__=='__main__':
    main()
