"""Converge a specified CuMnAs AC response; not the paper's intrinsic DC curve.

Run after pip install -e .: python examples/use_dynamic.py
Uses the original shifted [-pi,pi)^2 grid and a_nm=1 explicitly.
"""
from dataclasses import asdict
import json
from nmee_lab import DynamicSpec, evaluate_dynamic, cumnas_dynamic_model
from nmee_lab.ptbackend.cumnas import Model, points


def main():
    params = Model(a_nm=1.)
    model = cumnas_dynamic_model(theta_rad=.15, params=params)
    spec = DynamicSpec(.043+.017j, .071+.023j, .03, 100., 'A_z', 'x', 0)
    rows = []
    for n in (241, 481, 721):
        result = evaluate_dynamic(model, points(n)/params.a_nm, spec)
        rows.append(dict(n=n, response=result['response'], parts=result['parts']))
    residual = abs(rows[-1]['response']-rows[-2]['response'])/abs(rows[-1]['response'])
    output = dict(model='CuMnAs Eq5', theta_rad=.15, a_nm=params.a_nm,
                  settings=asdict(spec), convergence=rows,
                  protocol=result['protocol'], normalization=result['normalization'],
                  output_energy_eV=result['output_energy_eV'],
                  last_relative_mesh_change=residual,
                  scope='Full AC response at fixed switching rates; no eta-to-zero or intrinsic-DC claim')
    print(json.dumps(output, indent=2, default=lambda z: {'real':z.real,'imag':z.imag}, allow_nan=False))
    if residual > 1e-5:
        raise SystemExit('Mesh tolerance 1e-5 not reached for this benchmark')


if __name__ == '__main__':
    main()
