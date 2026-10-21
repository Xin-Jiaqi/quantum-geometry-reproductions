#!/usr/bin/env python3
"""v54: verified geometry from external H(k), analytic dH/dk NPZ or test lattice.

Input .npz keys hamiltonian=(...,N,N), velocities=(2,...,N,N);
no pickle and no assumed original-paper transport formula.
"""
import argparse,json
from pathlib import Path
import numpy as np
from qgmap.orbital_geometry_v54 import (
    band_geometry,periodic_three_band_model,moment_curvature_kernel,curvature_moment_product)

def evaluate(npz=None,N=36,mu=-1.1,temperature=.45):
    if npz is None:
        _,H,V=periodic_three_band_model(N)
        src={'type':'analytic_periodic_three_band_reference','N':N}
    else:
        with np.load(npz,allow_pickle=False) as obj:
            if not {'hamiltonian','velocities'} <= set(obj.files):
                raise ValueError('NPZ requires hamiltonian and velocities arrays')
            H=np.asarray(obj['hamiltonian'])
            V=np.asarray(obj['velocities'])
        src={'type':'user_H_and_dHdk_npz','path':str(npz)}
    g=band_geometry(H,V)
    kernel=moment_curvature_kernel(g,mu,temperature)
    split=curvature_moment_product(g)
    a=kernel['cross_m_ne_l'];b=kernel['full']
    return {'version':'v54','status':g.acceptance,'source':src,
       'band_count':int(H.shape[-1]),'mesh_shape':list(H.shape[:-2]),
       'gap_min':g.min_direct_gap,
       'berry_trace_max_abs':float(np.max(abs(np.sum(g.curvature,axis=-1)))),
       'metric_min_diag':float(np.min(np.diagonal(g.metric,axis1=-2,axis2=-1))),
       'moment_curvature_sector':kernel,
       'cross_to_full_magnitude':float(abs(a/b)) if b else None,
       'decomposition_status':split['acceptance'],
       'limitations':['Not the full orbital MNHE or PHE',
           'Requires matched analytic velocities and k-independent orbital basis',
           'No non-Abelian degenerate-band solver',
           'Missing orbital embedding, positional shift and full 3-band diagrams',
           'BZ mean requires uniformly weighted periodic k mesh; no SI units']}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-npz',type=Path)
    p.add_argument('--N',type=int,default=36)
    p.add_argument('--mu',type=float,default=-1.1)
    p.add_argument('--temperature',type=float,default=.45)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    txt=json.dumps(evaluate(args.input_npz,args.N,args.mu,args.temperature),
                   indent=2,allow_nan=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(txt,encoding='utf8')
    else:print(txt,end='')
if __name__=='__main__':main()
