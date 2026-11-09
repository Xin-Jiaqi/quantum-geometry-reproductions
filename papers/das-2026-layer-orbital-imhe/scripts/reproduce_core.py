#!/usr/bin/env python3
"""Recompute the declared periodic layer/orbital benchmark; no author data used."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from imhe_lab import (mesh, layer_fields, geometry, transport, exact_transport,
                      mixed_difference, projector_curvature, exact_geometry)

ROOT = Path(__file__).resolve().parents[1]


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {k: serial(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [serial(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def reproduce():
    k,w=mesh(32);fields=layer_fields(k);g=geometry(*fields.jets())
    r=transport(g,w,.2,.15)
    broken=layer_fields(k,magnetic_order=.13)
    rb=transport(geometry(*broken.jets()),w,.2,.15)
    k2,w2=mesh(48);r2=transport(geometry(*layer_fields(k2).jets()),w2,.2,.15)
    steps=[]
    for h in (.004,.002,.001):
        hi=mixed_difference(lambda a,b:exact_transport(fields,w,.2,.15,a,b)[0],h,h)
        dr=mixed_difference(lambda a,b:exact_transport(broken,w,.2,.15,a,b)[1],h,h)
        steps.append(dict(step=h,intrinsic_finite_difference=float(hi),
                          intrinsic_abs_error=float(abs(hi-r['intrinsic_full'])),
                          drift_max_abs_error=float(np.max(abs(dr-rb['band_drift_full'])))))
    points=np.array([[.73,-.48],[1.21,.32]])
    fg=layer_fields(points)
    projector=projector_curvature(lambda q:layer_fields(q).at(.03,.04)[0],points)
    spectral=exact_geometry(*fg.at(.03,.04))[1]
    result=dict(
        status='MODEL_READY',
        protocol='fixed_reference_A5_moment_H0_plus_uP_minus_bM_periodic_3band',
        source='arXiv:2604.20249v1; Eq1, Eq6/Eq7 partial sectors, A16 product rule',
        units='e=hbar=1; dimensionless k; explicit torus d2k/(2pi)^2; per probe u*b; drift per tau; NOT calibrated SI',
        parameters=dict(background_gate=.23,mu=.2,thermal_energy=.15,mesh=32,
                        magnetic_order_for_nonzero_drift=.13),
        time_reversal_symmetric=r,
        broken_time_reversal=rb,
        validation=dict(field_steps=steps,
                        intrinsic_mesh32_to48_abs=float(abs(r2['intrinsic_full']-r['intrinsic_full'])),
                        projector_vs_spectral_max_abs=float(np.max(abs(projector-spectral))),
                        gapped_T0=transport(g,w,-.5,0)),
        source_code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted((ROOT/'src').rglob('*.py'))},
        limits=['Independent synthetic periodic model, not R5G or its 0.05 e^2/h peak',
                'Reference A5 moment frozen under probe fields; optional C declares extra mixed operator',
                'Full means stated Berry-current / band-drift functionals, not all orbital transport',
                'No Lorentz/phase-space drift/disorder vertices, arbitrary Wannier embedding or degeneracy solver',
                'No generic gate-odd/even assertion and no nonquantized insulating DC plateau imposed'])
    return serial(result),g,w


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path)
    args=parser.parse_args()
    result,g,w=reproduce()
    if args.outdir:
        args.outdir.mkdir(parents=True,exist_ok=True)
        (args.outdir/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        with (args.outdir/'mu_scan.csv').open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['mu','intrinsic_eq6','intrinsic_full_A16','sea','energy_surface'])
            for mu in np.linspace(-1.9,2.2,101):
                r=transport(g,w,float(mu),.15)
                writer.writerow([mu,r['intrinsic_eq6'],r['intrinsic_full'],
                                 r['intrinsic_terms']['sea'],r['intrinsic_terms']['energy_surface']])
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
