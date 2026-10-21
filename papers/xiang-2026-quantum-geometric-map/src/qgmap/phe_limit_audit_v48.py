"""v48: controlled original Eq.(11) PHE limit tests on *named* frozen-H0 RTA.

The source paper's Supplemental Material collision prescription has NOT been
retrieved/verified. This module is a quantitative non-equivalence-aware audit
of the original main-text symmetric spin PHE geometric function, not a claim
of author-method identity, physical Hall voltage, or exact Fig.2(b) units.

Models are evaluated at fixed mu,T,beta, UV while the *external* one-body RTA
eta tends to zero. An eta^2 extrapolation is a diagnostic; coefficients of the
subleading correction depend on this regulator and collision-target choice.
"""
from __future__ import annotations
from dataclasses import asdict
import numpy as np
from .dirac_mixed_rta_v46 import Spec
from .dirac_source_tensor_v47 import tensor, derivative_source_check


def _simple(spec:Spec,step:float):
    z=tensor(spec,step=step)
    return {'beta':spec.beta,'mu':spec.mu,'T':spec.temperature,'eta':spec.eta,
            'UV_cutoff':spec.cutoff,'eta_symmetric':z['eta_times_symmetric'],
            'eta_antisymmetric':z['eta_times_antisymmetric'],
            'paper_main_eq11_tanh':z['printed_eq11_geometry'],
            'ratio_to_main_eq11':z['ratio_to_main_eq11'],
            'xy':z['xy']['total'],'yx':z['yx']['total']}


def eta_limit(spec=Spec(beta=.125,mu=.3,temperature=.03,nr=140,nphi=160),
              etas=(.016,.008,.004,.002,.001),step=2e-5):
    if len(etas)<4 or not all(np.isfinite(x) and x>0 for x in etas):
        raise ValueError('require >=4 finite positive eta')
    if len(set(etas))!=len(etas) or not all(a>b for a,b in zip(etas,etas[1:])):
        raise ValueError('eta must be distinct and descending')
    if spec.beta<=0 or spec.mu==0 or step<=0:
        raise ValueError('nonzero doping/positive beta/step required for ratio')
    vals=[_simple(Spec(spec.beta,spec.mu,spec.temperature,float(e),
                      spec.cutoff,spec.nr,spec.nphi),step) for e in etas]
    xs=np.asarray(etas)**2
    ratios=np.array([z['ratio_to_main_eq11'] for z in vals])
    # Fitted only in the small eta region: n=3,4,5 points provide independent
    # truncation-sensitivity checks. Do not fit a free multiplicative factor.
    fits=[]
    for n in (3,4,5):
        if n>len(etas):continue
        slope,intercept=np.polyfit(xs[-n:],ratios[-n:],1)
        fits.append({'smallest_n':n,'eta_sq_slope':float(slope),
                     'ratio_intercept_eta0':float(intercept),
                     'absolute_intercept_minus_unity':float(abs(intercept-1))})
    anti_scaled=[float(z['eta_antisymmetric']/z['eta']**2) for z in vals]
    return {'protocol':'v48_frozen_H0_RTA_original_symmetric_PHE_eta2_extrapolation',
            'scientific_status':'PARTIAL','spec':asdict(spec),'values':vals,
            'fit_eta_squared':fits,'antisymmetric_eta_scaled_by_eta2':anti_scaled,
            'limitations':['The original author SM collision operator is not accessible and has NOT been reconstructed',
                           'Numerical extrapolation is conditional on this specified one-body RTA',
                           'Eq8 signed prefactor and absolute Fig2b ~2pi discrepancy remain OPEN']}


def chemical_potential_temperature(spec=Spec(beta=.125,mu=.3,temperature=.03,
        eta=.002,nr=140,nphi=160),mus=(.15,.3),temperatures=(.02,.055),step=2e-5):
    if not mus or not temperatures or spec.beta<=0 or step<=0:
        raise ValueError('nonempty mus/temps and positive beta/step required')
    out=[]
    for T in temperatures:
        for mu in mus:
            if not np.isfinite([mu,T]).all() or mu<=0 or T<=0:raise ValueError('positive mu/T')
            pos=_simple(Spec(spec.beta,mu,T,spec.eta,spec.cutoff,spec.nr,spec.nphi),step)
            neg=_simple(Spec(spec.beta,-mu,T,spec.eta,spec.cutoff,spec.nr,spec.nphi),step)
            out.append({'mu_abs':mu,'temperature':T,
                        'positive_mu':pos,'negative_mu':neg,
                        'odd_parity_residual':pos['eta_symmetric']+neg['eta_symmetric'],
                        'positive_ratio':pos['ratio_to_main_eq11'],
                        'negative_ratio':neg['ratio_to_main_eq11']})
    return {'protocol':'v48_frozen_H0_RTA_eq11_mu_temperature_parity',
        'scientific_status':'PARTIAL','results':out,
        'warning':'tanh finite-T extension is benchmark of the declared model, not proof of source supplemental RTA'}


def cutoff_sweep(spec=Spec(beta=.125,mu=.3,temperature=.03,eta=.002,
          nr=140,nphi=160),cutoffs=(1.,1.5,2.5,4.),step=2e-5):
    if len(cutoffs)<2 or not all(np.isfinite(x) and x>0 for x in cutoffs) or step<=0:
        raise ValueError('two positive UV cutoffs required')
    vals=[_simple(Spec(spec.beta,spec.mu,spec.temperature,spec.eta,float(uv),
                     spec.nr,spec.nphi),step) for uv in cutoffs]
    ratio=[x['ratio_to_main_eq11'] for x in vals]
    return {'protocol':'v48_frozen_H0_RTA_original_PHE_UV_stability',
        'scientific_status':'PARTIAL','values':vals,
        'max_ratio_spread':float(max(ratio)-min(ratio)),
        'warning':'Continuum radial cutoff test does not replace lattice UV completion or source-SM collision check'}


def full_audit():
    sp=Spec(beta=.125,mu=.3,temperature=.03,eta=.002,nr=140,nphi=160)
    return {'protocol':'v48_original_phe_robustness_and_source_provenance',
        'scientific_status':'PARTIAL',
        'source_evidence':{'main_text':'arXiv:2510.02661v2, main Eqs. 3-8,10-11 and closing RTA statement',
          'supplemental_material':'NOT independently accessible from APS or arXiv through current tool; no author collision-operator claims',
          'theoretical_scheme':'Researcher-defined frozen-H0 scalar RTA, corrected v47 electric source and symmetric tensor'},
        'eta_limit':eta_limit(sp),
        'mu_temperature':chemical_potential_temperature(sp),
        'cutoff':cutoff_sweep(sp),
        'finite_difference_Fermi_source':derivative_source_check(Spec(nr=18,nphi=32),2e-5),
        'overall_scientific_acceptance':'PARTIAL'}
