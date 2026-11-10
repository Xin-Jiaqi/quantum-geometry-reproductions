#!/usr/bin/env python3
"""Print source-audited TB Table S1 and local SHG divergent-term cancellation as JSON.

Run from this paper folder with PYTHONPATH=src, or via root scripts/core.py.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from deep_dbc import EPSILONS, SOURCE_E0, SOURCE_VRATIO, benchmark, nodal_parameters


def calculate():
    nodal = [nodal_parameters(ep) for ep in EPSILONS]
    cases = [benchmark(ep, q) for ep in (.75, .90, 1.03, 1.50, 2.0, 3.0)
             for q in (1e-3, 1e-4, 1e-5)]
    error_e0 = max(abs(x['epsilon0']-ref) for x,ref in zip(nodal,SOURCE_E0))
    error_v = max(abs(x['vv_over_vc']-ref) for x,ref in zip(nodal,SOURCE_VRATIO))
    return {
        'source': 'N. Zou et al., Physical Review Letters 135, 226901 (2025), DOI 10.1103/sqwy-xhdf',
        'protocol': 'SM (S6)+(S19), with leading C-singular from (S12); xxy, generic oblique approach',
        'observable': 'local SHG singularity budget, NOT full finite SHG or GeTe',
        'units': 'dimensionless k (lattice units), energy per SM TB; hbar=1; common e^3/(2*hbar^2) omitted',
        'epsilon1_nodal_parameters': nodal,
        'table_s1_max_absolute_epsilon0_error': error_e0,
        'table_s1_max_absolute_velocity_ratio_error': error_v,
        'shg_cancellation': cases,
        'validation': 'MODEL_READY: leading singular terms and 3-orbital TB only',
        'limits': [
            'No full SHG spectrum: finite S7 metric/curvature/velocity contributions are omitted',
            'Frequencies are kept off optical poles (omega=0.7*epsilon0); no disorder or dephasing',
            'No absolute muA/V^2 calibration, GeTe DFT/Wannier or SOC/HSE reproduction',
            'The singular C expression assumes one valence band and a near-degenerate conduction doublet',
        ],
    }


def export(data, directory):
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'validation.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    rows=data['shg_cancellation']
    with (directory/'shg_cancellation.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(8,3.2))
    groups=[(1.03,[x for x in rows if x['epsilon1']==1.03]),(.75,[x for x in rows if x['epsilon1']==.75])]
    for eps, subset in groups:
        qs=[x['q'] for x in subset]
        axes[0].loglog(qs,[x['inter_abs'] for x in subset],marker='o',lw=1,label=f'inter, eps1={eps}')
        axes[0].loglog(qs,[x['mixed_abs'] for x in subset],marker='x',lw=1,ls='--',label=f'mixed, eps1={eps}')
        axes[1].loglog(qs,[x['residual_abs'] for x in subset],marker='o',label=f'finite residual, eps1={eps}')
    axes[0].set(xlabel='Distance q to node (lattice units)',ylabel='|local SHG term| (reduced units)',title='Individual 1/q divergences')
    axes[1].set(xlabel='Distance q to node (lattice units)',ylabel='|inter + mixed singular|',title='Leading divergences cancel')
    for ax in axes: ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
    fig.tight_layout();fig.savefig(directory/'shg_cancellation.png',dpi=180);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path)
    args=parser.parse_args()
    data=calculate()
    if args.outdir:export(data,args.outdir)
    print(json.dumps(data,ensure_ascii=False,allow_nan=False))


if __name__=='__main__': main()

