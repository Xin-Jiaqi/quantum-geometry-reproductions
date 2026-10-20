#!/usr/bin/env python3
"""Round 6 independent d-vector / Gauss quadrature / TR-frequency audit.

Immutable historical checkpoints are not changed. This script writes *only*
round6_* raw/figure products. The figure/source graph is NOT used as fit input.
"""
from __future__ import annotations
import json,sys,csv,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from igmc_lab.d_vector_independent import (dvector_dirac,pauli_trace_geometry,
                                           midpoint_dirac,aggregate_response)
from igmc_lab.core import (dirac_hamiltonian,pair_kernels,dirac_response,
                           physical_response_si,energy_from_frequency)


def gauss_square(N,mu,photon_eV,*,alpha='x',spin='x',mass_sign=1):
    """Gauss-Legendre tensor product over full square BZ, no eigenvectors.

    For insulating mu=0 the integrand is analytic. Metallic mu is a
    discontinuous occupation domain and is NOT claimed converged here.
    """
    roots,weights=np.polynomial.legendre.leggauss(N)
    kx,ky=np.meshgrid(np.pi*roots,np.pi*roots,indexing='ij')
    wx,wy=np.meshgrid(weights,weights,indexing='ij')
    d,dgrad=dvector_dirac(kx,ky,mass_sign=mass_sign)
    q=pauli_trace_geometry(d,dgrad,alpha,spin)
    E=q['gap']/2
    delta=((mu>-E).astype(float)-(mu>E).astype(float))
    den=1-(photon_eV/q['gap'])**2
    # integral d²k/(2π)² = (1/4)*integral d²x on [-1,1]²
    weight=wx*wy/4
    C=np.sum(weight*delta*q['Z01'])
    D=np.sum(weight*delta*q['D_per_omega01'])
    K=np.sum(weight*delta*(q['Z01']+1j*photon_eV*q['D_per_omega01'])/den)
    return {'N':N,'mu_eV':mu,'C_Ainv':float(C),'D_peromega_Ainv_eVinv':float(D),
            'Kubo_real_Ainv':float(K.real),'Kubo_imag_Ainv':float(K.imag)}


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--no-write', action='store_true', help='check evidence without editing any stored file')
    opts=arg.parse_args()
    params=json.loads((ROOT/'parameters_round6.json').read_text())
    photon=energy_from_frequency(params['frequency_hz'],protocol=params['frequency_protocol'])
    assert params['g_factor']==2.0 and params['chemical_potential_eV']==0.0
    assert params['lattice_a_Angstrom']==1.0 and params['beta_eV']==0.25
    assert params['vF_eV_Angstrom']==0.825 and params['M_eV_Angstrom_squared']==0.1
    rng=np.random.default_rng(20862025)
    ka=rng.uniform(-np.pi,np.pi,(23,19));kb=rng.uniform(-np.pi,np.pi,(23,19))
    max_errors={'Z_abs':0.,'Dperomega_abs':0.}
    for sign in (1,-1):
        H,v=dirac_hamiltonian(ka,kb,mass_sign=sign)
        d,grad=dvector_dirac(ka,kb,mass_sign=sign)
        for alpha in ('x','y'):
            for spin in ('x','y','z'):
                q=pair_kernels(H,v,alpha,spin)
                analytic=pauli_trace_geometry(d,grad,alpha,spin)
                max_errors['Z_abs']=max(max_errors['Z_abs'],float(np.max(abs(q['Z01']-analytic['Z01']))))
                max_errors['Dperomega_abs']=max(max_errors['Dperomega_abs'],float(np.max(abs(2*q['Q01']/q['gap']-analytic['D_per_omega01']))))
    if max(max_errors.values())>1e-12:raise AssertionError(f'independent local identity failed: {max_errors}')
    # Full complex K(omega) and TR-related magnetic partners; each of the
    # six tensor components computed for two masses and both frequencies.
    complex_tests=[];max_tr=0.;max_spec=0.
    for frequency_hz in (1e12,3e12):
        omega=energy_from_frequency(frequency_hz,protocol='cyclic_hz')
        for alpha in ('x','y'):
            for spin in ('x','y','z'):
                x=midpoint_dirac(81,0.,omega,alpha=alpha,spin=spin,mass_sign=1)
                y=midpoint_dirac(81,0.,omega,alpha=alpha,spin=spin,mass_sign=-1)
                # TR conjugates the analytic complex response at positive omega:
                # C[TR(H)] = C[H], D[TR(H)] = -D[H].
                tr=max(abs(x['Kubo_real']-y['Kubo_real']),abs(x['Kubo_imag']+y['Kubo_imag']))
                max_tr=max(max_tr,float(tr))
                ref=dirac_response(81,0.,omega_eV=omega,alpha=alpha,beta=spin)
                spectral=max(abs(x['Kubo_real']-ref['Kubo_real']),abs(x['Kubo_imag']-ref['Kubo_imag']))
                max_spec=max(max_spec,float(spectral))
                complex_tests.append({'frequency_hz':frequency_hz,'alpha':alpha,'spin':spin,
                    'C_real':x['Kubo_real'],'D_quadrature':x['Kubo_imag'],
                    'tr_C_real':y['Kubo_real'],'tr_D_quadrature':y['Kubo_imag'],
                    'TR_conjugation_error':float(tr),'independent_vs_spectral_error':float(spectral)})
    if max_tr>1e-10 or max_spec>1e-10:raise AssertionError('TR or full Kubo response relation failed')

    # This is an independent *integration scheme* from all older midpoint
    # grids, and includes the full frequency dependence without matrix solvers.
    gauss=[gauss_square(n,params['chemical_potential_eV'],photon) for n in params['gauss_legendre_orders']]
    midpoint=midpoint_dirac(params['midpoint_reference_N'],0.,photon,spin='x')
    src=dirac_response(params['midpoint_reference_N'],0.,omega_eV=photon,alpha='x',beta='x')
    dist=abs(gauss[-1]['D_peromega_Ainv_eVinv']-gauss[-2]['D_peromega_Ainv_eVinv'])
    if dist>2e-8:raise AssertionError(f'Gauss integration not converged: {dist}')
    if abs(gauss[-1]['D_peromega_Ainv_eVinv']-midpoint['sigma_D_per_omega'])>2e-6:
        raise AssertionError('independent Gauss and midpoint differ')
    si_main=physical_response_si(gauss[-1]['D_peromega_Ainv_eVinv']*photon,g_factor=params['g_factor'])
    si_twice=2*si_main
    record={'parameters_round6':params,'protocol':'PRL Eq.(17) d-vector identity -> no eigenvectors -> independent Gauss quadrature',
        'source_fig2d_values':'manual approximate only; no author raw data; neither source curve used as input',
        'photon_eV':photon,'pointwise_checks':12,'max_pair_kernels_abs_error':max_errors,
        'TR_exact_complex_checks':len(complex_tests),'max_complex_TR_conjugation_error':max_tr,
        'max_complex_vs_spectral_Kubo_error':max_spec,
        'gauss_convergence_abs_Dperomega_last_step':dist,
        'gauss_orders':[x['N'] for x in gauss],
        'gauss_Dperomega': [x['D_peromega_Ainv_eVinv'] for x in gauss],
        'midpoint_N301_Dperomega':midpoint['sigma_D_per_omega'],
        'spectral_N301_Dperomega':src['sigma_D_per_omega'],
        'gauss_main_fig2d_gap_SI_A_m1_T1':si_main,
        'SM_eq34_doubled_intermediate_diagnostic_only':si_twice,
        'conclusion':'Independent 3-vector formulas and Gauss quadrature support main Eq.(12). Author plot absolute normalization remains unresolved.'}
    if not opts.no_write:
        raw=ROOT/'results/raw';raw.mkdir(exist_ok=True,parents=True)
        (raw/'round6_analytic_dirac_audit.json').write_text(json.dumps(record,indent=2)+'\n')
        with (raw/'round6_gauss_convergence.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=gauss[0].keys());w.writeheader();w.writerows(gauss)
        with (raw/'round6_complex_TR_quadrature.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=complex_tests[0].keys());w.writeheader();w.writerows(complex_tests)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ns=np.array([p['N'] for p in gauss])
        errors=np.array([abs(p['D_peromega_Ainv_eVinv']-gauss[-1]['D_peromega_Ainv_eVinv']) for p in gauss[:-1]])
        fig,ax=plt.subplots(figsize=(7.4,4.6),layout='constrained')
        ax.semilogy(ns[:-1],errors,'o-',label='Gauss-Legendre vs N=128 reference')
        ax.axhline(abs(midpoint['sigma_D_per_omega']-gauss[-1]['D_peromega_Ainv_eVinv']),
                   linestyle='--',label='Independent midpoint N=301 deviation')
        ax.set(xlabel='Gauss-Legendre order along each BZ direction',
               ylabel='Absolute error in D/omega [1/(Angstrom eV)]',
               title='Dirac in-gap displacement: independent BZ quadrature')
        ax.text(.03,.08,f'1 THz, beta=0.25 eV, g=2\nGauss result: {si_main:.8f} A/(m T)\nNo figure-derived rescaling',
                transform=ax.transAxes,fontsize=9,
                bbox=dict(facecolor='white',alpha=.9,edgecolor='.8'))
        ax.legend(fontsize=9);ax.grid(alpha=.22)
        figdir=ROOT/'results/figures';figdir.mkdir(exist_ok=True,parents=True)
        fig.savefig(figdir/'round6_independent_gauss_convergence.png',dpi=170)
        plt.close(fig)
    print(json.dumps(record,indent=2))
    return record

if __name__=='__main__':main()
