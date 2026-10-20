#!/usr/bin/env python3
"""Round 7: Fermi-contour metal audit + independent Eq.(13) + release matrix.

Runs the paper-printed Eq.(11)/(12) branch with no broadening, no fitted
normalization. Does not write when --no-write is supplied.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from igmc_lab.core import (energy_from_frequency,physical_response_si)
from igmc_lab.fermi_contour import (integrate_contour,integrate_conduction_fermi_line,contour_roots_x)

def main(argv=None):
    ap=argparse.ArgumentParser();ap.add_argument('--no-write',action='store_true');x=ap.parse_args(argv)
    pin=json.loads((ROOT/'parameters_round7.json').read_text())
    photon=energy_from_frequency(pin['frequency_hz'],protocol=pin['frequency_protocol'])
    assert pin['model'] == 'paper_dirac_eq17' and pin['temperature_K']==0
    mu_list=pin['convergence_mu_eV']
    rows=[]; max_c_err=0.;max_d_err=0.;max_fs_err=0.;max_mirror_err=0.
    for mu in mu_list:
        outputs={}
        for n in pin['orders']:
            z=integrate_contour(mu,outer_order=n,inner_order=pin['inner_order'],photon_eV=photon)
            outputs[n]=z
        fs=integrate_conduction_fermi_line(mu,order=pin['fs_order'])
        a,b=outputs[pin['orders'][-2]],outputs[pin['orders'][-1]]
        max_c_err=max(max_c_err,abs(b['C_xy']-a['C_xy']))
        max_d_err=max(max_d_err,abs(b['D_peromega_xx']-a['D_peromega_xx']))
        max_fs_err=max(max_fs_err,abs(b['C_xy']-fs['C_xy_eq13']))
        neg=integrate_contour(-mu,outer_order=pin['orders'][-1],inner_order=pin['inner_order'],photon_eV=photon)
        max_mirror_err=max(max_mirror_err,abs(b['C_xy']-neg['C_xy']),abs(b['D_peromega_xx']-neg['D_peromega_xx']))
        for n in pin['orders']:
            z=outputs[n]
            rows.append({'mu_eV':mu,'outer_order':n,'inner_order':pin['inner_order'],
                 'C_xy_Ainv':z['C_xy'],'D_peromega_xx_Ainv_eVinv':z['D_peromega_xx'],
                 'C_xy_SI_A_per_mT':physical_response_si(z['C_xy']),
                 'D_xx_SI_A_per_mT':physical_response_si(z['D_peromega_xx']*photon),
                 'Kubo_C_xy_SI':physical_response_si(z['K_real_xy']),
                 'Kubo_D_xx_SI':physical_response_si(z['K_imag_xx']),
                 'independent_fermi_line_C_xy_Ainv':fs['C_xy_eq13'],
                 'fs_minus_C_Ainv':fs['C_xy_eq13']-z['C_xy'],
                 'ky_separatrix_count':len(z['topology_splits_ky'])-2,
                 'max_x_FS_crossings':z['max_x_crossings']})
    # Finite-frequency low-omega errors; do not compare to author's plotted sources.
    main_mu=1.;part=integrate_contour(main_mu,outer_order=pin['orders'][-1],
        inner_order=pin['inner_order'],photon_eV=photon)
    freq_relative_C=abs(part['K_real_xy']/part['C_xy']-1)
    freq_relative_D=abs(part['K_imag_xx']/(part['D_peromega_xx']*photon)-1)
    # Curve shows location of peak and electron/hole symmetry without using publisher pixels.
    curve=[]
    for mu in pin['scan_mu_eV']:
        z=integrate_contour(mu,outer_order=pin['curve_order'],inner_order=pin['inner_order'],photon_eV=photon)
        curve.append({'abs_mu_eV':mu,'C_xy_Ainv':z['C_xy'],'D_peromega_xx_Ainv_eVinv':z['D_peromega_xx'],
                      'C_xy_SI_A_per_mT':physical_response_si(z['C_xy']),
                      'D_xx_SI_A_per_mT':physical_response_si(z['D_peromega_xx']*photon)})
    cpeak=max(curve,key=lambda r:abs(r['C_xy_SI_A_per_mT']))
    dpeak=max(curve,key=lambda r:abs(r['D_xx_SI_A_per_mT']))
    summary={'protocol':'Xi2025 PRL Eq11/12, T=0 Fermi domain partition, no eta',
      'parameters_file':'parameters_round7.json','photon_eV':photon,
      'tested_mu_eV':mu_list,'orders':pin['orders'],'fs_order':pin['fs_order'],
      'max_absolute_convergence_C_Ainv':max_c_err,
      'max_absolute_convergence_Dperomega_Ainv_eVinv':max_d_err,
      'max_absolute_occupied_vs_FS_eq13_C_Ainv':max_fs_err,
      'max_absolute_particle_hole_mirror_error':max_mirror_err,
      'mu1_C_xy_Ainv':part['C_xy'], 'mu1_Dperomega_xx_Ainv_eVinv':part['D_peromega_xx'],
      'mu1_C_SI_A_per_mT':physical_response_si(part['C_xy']),
      'mu1_D_SI_A_per_mT':physical_response_si(part['D_peromega_xx']*photon),
      'mu1_full_Kubo_relative_deviation_C':freq_relative_C,
      'mu1_full_Kubo_relative_deviation_D':freq_relative_D,
      'curve_C_peak':cpeak,'curve_D_peak':dpeak,
      'independent_checks':'Fermi-contour root partition and Eq13 1D Fermi-line delta integral',
      'figure_comparison':'visual estimate only, paper raw data not available; not used for fitting',
      'overall_paper_status':'PARTIAL'}
    assert len(contour_roots_x(.5,1.))==2
    if max_c_err>4.e-8 or max_d_err>2.e-8:
        raise AssertionError('not yet metal convergence: '+str((max_c_err,max_d_err)))
    if max_fs_err>2.e-8:raise AssertionError('occupied Fermi sea vs Eq13 contour mismatch')
    if max_mirror_err>1e-13:raise AssertionError('particle-hole symmetry mismatch')
    if max(freq_relative_C,freq_relative_D)>2e-4:
        raise AssertionError('frequency near-DC expansion not consistent')
    if not x.no_write:
        raw=ROOT/'results/raw';raw.mkdir(parents=True,exist_ok=True)
        with (raw/'round7_metal_convergence.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
        with (raw/'round7_metal_curve.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=curve[0].keys());w.writeheader();w.writerows(curve)
        (raw/'round7_metal_audit.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        out=ROOT/'results/figures';out.mkdir(parents=True,exist_ok=True)
        xs=np.array([r['abs_mu_eV'] for r in curve])
        for key,label,fname in [
          ('C_xy_SI_A_per_mT',r'Conduction: $-g\mu_B\sigma^C_{xy}$', 'round7_fig2b_metal_Fermi_boundary.png'),
          ('D_xx_SI_A_per_mT',r'Displacement: $-g\mu_B\sigma^D_{xx}$', 'round7_fig2d_metal_Fermi_boundary.png')]:
            fig,ax=plt.subplots(figsize=(7.5,4.8),layout='constrained')
            ys=np.array([r[key] for r in curve]);ax.plot(xs,ys,'o-',ms=2.6,lw=1.4,label='Exact Fermi boundary + Gauss quadrature')
            ax.set(xlabel=r'$|\mu|$ (eV); positive/negative identical in Eq. (17)',
                   ylabel='Current / magnetic field (A / (m T))',title=label)
            ax.axvline(1.,ls='--',lw=1.,c='.5',label='1 eV comparison point')
            ax.grid(alpha=.25);ax.legend(fontsize=8)
            ax.text(.03,.05,'Paper main Eqs. (11)-(12), 1 THz, g=2\nNo source-image amplitude fitting',
                    transform=ax.transAxes,fontsize=9,va='bottom')
            fig.savefig(out/fname,dpi=175);plt.close(fig)
        fig,ax=plt.subplots(figsize=(7.4,4.7),layout='constrained')
        selected=[r for r in rows if r['mu_eV']==1.0]
        ns=np.array([r['outer_order'] for r in selected]);ref=selected[-1]
        ax.semilogy(ns[:-1],[abs(r['C_xy_Ainv']-ref['C_xy_Ainv']) for r in selected[:-1]],'o-',label='C sea vs highest order')
        ax.semilogy(ns[:-1],[abs(r['D_peromega_xx_Ainv_eVinv']-ref['D_peromega_xx_Ainv_eVinv']) for r in selected[:-1]],'s-',label='D/omega sea vs highest order')
        ax.axhline(abs(ref['fs_minus_C_Ainv']),ls=':',label='Eq.(13) Fermi-line vs sea (highest order)')
        ax.set(xlabel='Split outer Gauss order (each Fermi-topology ky segment)',
               ylabel='Absolute dimensionful kernel difference',title='Metallic mu=1 eV: exact Fermi-domain convergence')
        ax.grid(alpha=.25);ax.legend(fontsize=8);fig.savefig(out/'round7_metal_convergence.png',dpi=175);plt.close(fig)
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    return summary
if __name__=='__main__':main()
