"""Find signed insulating zeros of bilayer KM Cxx under scalar and surface perturbations."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from run_fixed_stack_layer_cxx import (AuditBilayer, AuditSpec, ab_ba_shift, fixed_stack_matrices,
                                       integration_grid, response)


def study(args):
    model = AuditBilayer(AuditSpec(shift=ab_ba_shift(args.fraction), staggered_eV=.045, tunneling_eV=1.))
    k,w,_ = integration_grid(model, *args.grid)
    blocks = [fixed_stack_matrices(model, k[i:i+512]) for i in range(0,len(k),512)]
    h = np.concatenate([b[0] for b in blocks])
    v = np.concatenate([b[1] for b in blocks])
    del blocks
    # A local, nonmagnetic crystal field on the support-facing layer. This term is k-independent.
    sublattice = np.array([1,1,-1,-1])
    correction = np.concatenate((args.top_surface*sublattice,args.surface*sublattice))
    h += np.diag(correction) / 1000
    rows = {}

    def calc(u):
        if u not in rows:
            rows[u] = response(model,h,v,w,10.,u+args.offset,.03125)
        return rows[u]['total']

    grid = np.arange(args.low, args.high+.01*args.spacing, args.spacing)
    values = [calc(u) for u in grid]
    roots = []
    for u0,u1,y0,y1 in zip(grid,grid[1:],values,values[1:]):
        if y0 is None or y1 is None or y0*y1 > 0:
            continue
        def insulating_value(u):
            y = calc(u)
            if y is None:
                raise ValueError("metal in bracket")
            return y
        try:
            root = brentq(insulating_value,u0,u1,xtol=.001)
        except ValueError:
            continue
        if not roots or abs(root-roots[-1]['Uext'])>.03:
            roots.append(dict(Uext=root, bracket=[u0,u1], **rows[root]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data = dict(domain_fraction=args.fraction, surface_meV=args.surface, top_surface_meV=args.top_surface, offset_meV=args.offset,
                grid=args.grid, window=[args.low,args.high], spacing=args.spacing, roots=roots)
    Path(str(args.output)+'.json').write_text(json.dumps(data,indent=2)+'\n')
    with Path(str(args.output)+'.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=['U_ext_meV',*next(iter(rows.values()))],lineterminator='\n')
        writer.writeheader()
        for u,r in sorted(rows.items()):
            writer.writerow(dict(U_ext_meV=u,**r))
    print(args.output.name,[(round(r['Uext'],3),round(r['gap_meV'],2)) for r in roots],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fraction',type=float,required=True)
    p.add_argument('--surface',type=float,default=0.)
    p.add_argument('--top-surface',type=float,default=0.)
    p.add_argument('--offset',type=float,default=0.)
    p.add_argument('--grid',type=int,nargs=3,default=[31,30,48])
    p.add_argument('--low',type=float,default=-16.)
    p.add_argument('--high',type=float,default=16.)
    p.add_argument('--spacing',type=float,default=2.)
    p.add_argument('--output',type=Path,required=True)
    study(p.parse_args())
