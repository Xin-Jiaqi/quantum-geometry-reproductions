#!/usr/bin/env python3
"""v40 RA/AR/RR/AA QWZ noncrossing current-ladder diagnostic, NOT physical Hall."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

from qgmap.current_vertex_v40 import current_audit, exact_disorder_ra_torus


def evaluate(mesh=21,energy=1.1,W=.12,eta=.12,L=4,seeds=6):
    return {
        'schema':'qgmap-v40-model-1',
        'science_status':'PARTIAL',
        'validation_scope':'separate SCBA ladder channels + finite-torus bare RA',
        'SCBA':current_audit(mesh=mesh,energy=energy,W=W,eta=eta),
        'finite_torus':exact_disorder_ra_torus(L=L,W=W,mu=energy,eta=eta,seeds=seeds),
        'physical_caveat':'The antisymmetric RA correlator is not Hall conductivity; the 4-channel symmetric spectral sum cannot replace Kubo-Streda/Bastin Fermi-sea and surface terms.'
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',type=int,default=21)
    ap.add_argument('--energy',type=float,default=1.1)
    ap.add_argument('--W',type=float,default=.12)
    ap.add_argument('--eta',type=float,default=.12)
    ap.add_argument('--L',type=int,default=4)
    ap.add_argument('--seeds',type=int,default=6)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    data=json.dumps(evaluate(args.mesh,args.energy,args.W,args.eta,args.L,args.seeds),indent=2,allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(data+'\n',encoding='utf-8')
    else:print(data)
if __name__=='__main__':main()
