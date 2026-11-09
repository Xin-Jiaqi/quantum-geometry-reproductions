#!/usr/bin/env python3
"""Recompute Ezawa 2025 MIC benchmark, optionally exporting figures and CSVs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from altermagnet_bpve import (
    Params, analytic_xxx_eq17, analytic_xyy_printed_eq19,
    band_observables, contour_injection, grid_injection,
)


def calculate():
    p = Params()
    energies = np.linspace(.205, 1.80, 96)
    sharp = np.array([contour_injection(e, p, ntheta=720) for e in energies])
    approx = np.array([analytic_xxx_eq17(e, p) for e in energies])
    neel = np.linspace(-.18, +.18, 73)
    scans = np.array([contour_injection(.4, Params(altermagnet=float(j)), kmax=2., ntheta=720)
                      for j in neel])
    psmall = Params(altermagnet=.005)
    audit_xyy = contour_injection(.4, psmall, optical='yy', ntheta=720)
    printed_xyy = analytic_xyy_printed_eq19(.4, psmall)
    grid_xxx = grid_injection(.4, Params(altermagnet=.03), nk=651, kmax=1.2, gamma=.01)
    contour_xxx = contour_injection(.4, Params(altermagnet=.03))
    return {
        "paper": "Motohiko Ezawa, PRB 111, L201405 (2025), DOI 10.1103/PhysRevB.111.L201405",
        "protocol": "mic.sharp_contour.quantum_metric",
        "formula": "source Eq.(7)/(9)/(11) vs Eq.(17); symmetric LPL optical indices",
        "normalization": "sigma / [tau*e^3/(hbar^2*W*k0)], with positive e; hbar=e=tau=W*k0=1",
        "parameters": vars(p),
        "photon_energies": energies.tolist(),
        "sigma_x_xx_contour": sharp.tolist(),
        "sigma_x_xx_eq17": approx.tolist(),
        "max_relative_eq17_difference": float(np.max(np.abs((sharp-approx)/np.maximum(np.abs(approx), 1e-15)))),
        "J_scan": neel.tolist(),
        "sigma_x_xx_J_scan": scans.tolist(),
        "neel_reversal_max_odd_error": float(np.max(np.abs(scans+scans[::-1]))),
        "grid_vs_contour": {
            "J": .03, "photon_energy": .4, "nk": 651, "kmax": 1.2, "gamma": .01,
            "broadened": float(grid_xxx), "sharp": float(contour_xxx),
            "relative_difference": float(abs((grid_xxx-contour_xxx)/contour_xxx)),
        },
        "source_eq19_audit": {
            "component": "x;yy", "J": .005, "photon_energy": .4,
            "numeric_contour": float(audit_xyy),
            "literal_printed_eq19": float(printed_xyy),
            "relative_difference": float(abs((audit_xyy-printed_xyy)/printed_xyy)),
            "status": "OPEN_SOURCE_DISCREPANCY; no corrective factor fitted",
        },
        "scope": ["two-band continuum only", "single monotonic optical contour below Pauli cutoff",
                  "T=0 Fermi occupation mu=0", "model units; no SI/Wannier material integration"],
    }


def _write_csv(path: Path, names: list[str], columns):
    with path.open('w', newline='', encoding='utf8') as f:
        writer = csv.writer(f)
        writer.writerow(names)
        writer.writerows(zip(*columns))


def export_results(data, out: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    out.mkdir(exist_ok=True, parents=True)
    (out/'validation.json').write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')
    E = np.array(data['photon_energies']); exact = np.array(data['sigma_x_xx_contour'])
    analytical = np.array(data['sigma_x_xx_eq17'])
    J = np.array(data['J_scan']); values = np.array(data['sigma_x_xx_J_scan'])
    _write_csv(out/'mic_spectrum.csv', ['photon_energy_over_eps0', 'contour', 'eq17'], [E, exact, analytical])
    _write_csv(out/'mic_J_reversal.csv', ['J_over_eps0_k0minus2', 'sigma_x_xx'], [J, values])

    plt.rcParams.update({'font.size': 9, 'axes.linewidth': .75, 'savefig.dpi': 250})
    blue, ochre = '#236b95', '#ca892f'
    p = Params(); k=np.linspace(-2., 2., 600)
    fig, ax=plt.subplots(figsize=(5.2, 3.0))
    for j,color,label in [(0, '#85929e', 'J=0'), (.1, blue, 'J=0.1')]:
        s=band_observables(k,np.zeros_like(k),Params(altermagnet=j))
        ax.plot(k,s['upper'],color=color,lw=1.35,label=label)
        ax.plot(k,s['lower'],color=color,lw=1.35)
    ax.set(xlabel=r'$k_x/k_0$  ($k_y=0$)',ylabel=r'$E/\varepsilon_0$',title='Rashba + d-wave altermagnetic asymmetry')
    ax.legend(frameon=False,fontsize=8);ax.axvline(0,color='0.7',lw=.5)
    fig.tight_layout();fig.savefig(out/'fig_band_asymmetry.png');fig.savefig(out/'fig_band_asymmetry.svg');plt.close(fig)

    fig,ax=plt.subplots(figsize=(4,3.7));q=np.linspace(-.3,.3,501)
    x,y=np.meshgrid(q,q,indexing='ij')
    for j,color,label in [(0., '#85929e','J=0'),(.1,blue,'J=+0.1'),(-.1,ochre,'J=-0.1')]:
        spec=band_observables(x,y,Params(altermagnet=j))['gap']
        ax.contour(x,y,spec,levels=[.4],colors=[color],linewidths=[1.3])
        ax.plot([],[],color=color,label=label)
    ax.legend(frameon=False,fontsize=8)
    ax.set(xlabel=r'$k_x/k_0$',ylabel=r'$k_y/k_0$',title=r'Resonant contours $\Delta\varepsilon=0.4\varepsilon_0$')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(out/'fig_resonant_contour.png');fig.savefig(out/'fig_resonant_contour.svg');plt.close(fig)

    fig,ax=plt.subplots(figsize=(5.3,3.2))
    scale = -p.altermagnet/(8*p.rashba)
    ax.plot(E,exact/scale,color=blue,lw=1.6,label='Numerical contour')
    ax.plot(E,analytical/scale,color=ochre,lw=1,ls='--',label='PRB Eq. (17)')
    ax.axhline(1,color='0.6',ls=':',lw=.85,label='High-energy plateau')
    ax.set(xlabel=r'$\hbar\omega/\varepsilon_0$',ylabel=r'$\sigma^{x;xx}/[-\tau e^3 J/(8\hbar^2Wk_0\lambda)]$',title='Magnetic injection current (LPL)')
    ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/'fig_mic_spectrum.png');fig.savefig(out/'fig_mic_spectrum.svg');plt.close(fig)

    fig,ax=plt.subplots(figsize=(5.2,3.2))
    ax.plot(J,values,color=blue,lw=1.6,label='Numerical')
    ax.plot(J,[analytic_xxx_eq17(.4,Params(altermagnet=float(j))) for j in J],
            color=ochre,lw=1,ls='--',label='Eq. (17)')
    ax.axhline(0,color='0.7',lw=.7)
    ax.set(xlabel=r'$J/(\varepsilon_0/k_0^2)$',ylabel=r'$\sigma^{x;xx}/[\tau e^3/(\hbar^2Wk_0)]$',
           title='Néel-order reversal')
    ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/'fig_neel_reversal.png');fig.savefig(out/'fig_neel_reversal.svg');plt.close(fig)

    fig,ax=plt.subplots(figsize=(4.7,3.8))
    q=np.linspace(-.45,.45,501);X,Y=np.meshgrid(q,q,indexing='ij');ob=band_observables(X,Y,p)
    delta=.015/(np.pi*((ob['gap']-.4)**2+.015**2))
    kernel=ob['velocity_difference'][...,0]*ob['metric'][...,0,0]*delta
    im=ax.imshow(kernel.T,origin='lower',extent=[q.min(),q.max(),q.min(),q.max()],
                 cmap='RdBu_r',vmin=-np.max(abs(kernel)),vmax=np.max(abs(kernel)),aspect='equal')
    ax.set(xlabel=r'$k_x/k_0$',ylabel=r'$k_y/k_0$',title=r'$\Delta v^x\,g_{xx}\,\delta_\gamma(2d-0.4)$')
    fig.colorbar(im,ax=ax,shrink=.8,label='Model units')
    fig.tight_layout();fig.savefig(out/'fig_k_resolved_kernel.png');fig.savefig(out/'fig_k_resolved_kernel.svg');plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',type=Path, help='Write figures, CSVs, validation.json')
    args=ap.parse_args()
    results=calculate()
    if args.outdir:
        export_results(results,args.outdir)
    print(json.dumps(results,ensure_ascii=False,allow_nan=False))


if __name__=='__main__':
    main()
