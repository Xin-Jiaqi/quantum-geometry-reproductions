"""Stage 4: orbital-hopping interference + mixed Berry-curvature slope.

Reproduce with OPENBLAS_NUM_THREADS=1, installed local PRL reproduction.
Only 0K gapped static frozen-occupation response, not the PRL driven DC.
"""
from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
import argparse,json
import numpy as np
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.response_mechanism import evaluate_mechanism
from scripts.run_stage3_sign_reversal import spec_at
from scripts.sliding_static_scan import response_weighted


def run(out_path,quick=False):
    default=spec_at(0)
    G=SlidingBilayer(default).model().reciprocal
    k,w,meta=valley_quadrature(G,41,1.5,35,64)
    k2,w2,meta2=valley_quadrature(G,61,1.5,60,80)
    variations={
      'baseline':{},
      'weak_hopping_10meV':{'tunneling_eV':.010},
      'intermediate_hopping_15meV':{'tunneling_eV':.015},
      'weaker_rashba_15meV':{'rashba_soc_eV':.015},
      'stronger_intrinsic_SOC_12meV':{'intrinsic_soc_eV':.012},
      'stronger_intrinsic_SOC_14meV':{'intrinsic_soc_eV':.014},
      'longer_hop_decay_0p18nm':{'decay_nm':.18},
      'longer_hop_decay_0p20nm':{'decay_nm':.20},
    }
    if quick: variations={k:v for k,v in variations.items() if k in ('baseline','weak_hopping_10meV')}
    results={}
    for name,changes in variations.items():
        cache={}
        def ev(u):
            key=round(float(u),9)
            if key not in cache:cache[key]=evaluate_mechanism(replace(spec_at(u),**changes),k,w)
            return cache[key]
        left,right=.10,.40
        a,b=ev(left),ev(right)
        if not (a['gamma'] < 0 < b['gamma'] and a['sampled_min_mu_distance_eV']>.0002 and b['sampled_min_mu_distance_eV']>.0002):
            raise AssertionError(f'Not a gapped sign reversal in perturbation {name}: {a}, {b}')
        vals={str(x):ev(x) for x in (.1,.2,.4)}
        for _ in range(7 if quick else 12):
            mid=(left+right)/2
            e=ev(mid)
            if e['gamma']>0:right=mid
            else:left=mid
        results[name]={'parameter_changes':changes,
                       'endpoints_and_midpoint':vals,
                       'zero_bracket':[left,right],
                       'zero_midpoint':(left+right)/2,
                       'min_sampled_mu_distance_eV':min(x['sampled_min_mu_distance_eV'] for x in cache.values())}
        print('robust',name,'root=',(left+right)/2,'mu_min=',results[name]['min_sampled_mu_distance_eV'],flush=True)
    # Orbital-channel inclusion-exclusion measures NONADDITIVITY; not a unique band-gauge contribution.
    channels={}
    for u in (.05,.1,.15,.2,.3,.4,1.0):
        s=spec_at(u)
        c={channel:evaluate_mechanism(s,k,w,hopping_channel=channel)
           for channel in ('full','same_only','cross_only')}
        c['inclusion_exclusion_nonadditivity']=c['full']['gamma']-c['same_only']['gamma']-c['cross_only']['gamma']
        channels[str(u)]=c
    # Independently validate spectral expression against original numerical kernel.
    benchmarks={}
    for u in (.1,.4):
        s=spec_at(u)
        old=response_weighted(s,k,w)['response']['total_x']
        new=channels[str(u)]['full']['gamma']
        if abs(old-new)>1e-12:raise AssertionError('new Berry-kernel does not match frozen original')
        fine=evaluate_mechanism(s,k2,w2)
        if abs(fine['gamma']-new)>5e-9:raise AssertionError('mechanism diagnostic did not converge')
        benchmarks[str(u)]={'prior_independent_kernel':old,'stage4_full':new,
                             'high_resolution':fine,'mesh_change':abs(fine['gamma']-new)}
    # Relative hopping phase: a new physical Hamiltonian, NOT a basis change.
    relative_channel_phase={}
    for u in (.1,.2,.4):
        s=spec_at(u)
        normal=evaluate_mechanism(s,k,w)['gamma']
        reversed_cross=evaluate_mechanism(s,k,w,cross_amplitude_sign=-1)['gamma']
        relative_channel_phase[str(u)]={'cross_plus_gamma':normal,
                                      'cross_minus_gamma':reversed_cross,
                                      'odd_under_cross_sign':(normal-reversed_cross)/2,
                                      'even_under_cross_sign':(normal+reversed_cross)/2}
    # Layer-gauge sign on all interlayer hoppings leaves every layer-diagonal spin response invariant.
    layer_gauge_evenness={}
    for t in (.0025,.005,.01,.02):
        s=replace(spec_at(.1),tunneling_eV=t)
        pos=evaluate_mechanism(s,k,w)['gamma']
        neg=evaluate_mechanism(replace(s,tunneling_eV=-t),k,w)['gamma']
        if abs(pos-neg)>5e-12:raise AssertionError('layer-gauge parity violated')
        layer_gauge_evenness[str(t)]={'positive_t':pos,'negative_t':neg,'even_residual':pos-neg,'gamma_over_t_sq':pos/t**2}
    # Same vs cross in weak hopping: genuinely non-additive, not just interpolation at full strength.
    weak={}
    for t in (.005,.01,.015,.020):
        q={}
        for u in (.1,.4):
            vals={ch:evaluate_mechanism(replace(spec_at(u),tunneling_eV=t),k,w,hopping_channel=ch)['gamma']
                  for ch in ('full','same_only','cross_only')}
            vals['inclusion_exclusion_nonadditivity']=vals['full']-vals['same_only']-vals['cross_only']
            q[str(u)]=vals
        weak[str(t)]=q
    root=ROOT/'results'
    report={
      'model':'eight-band bilayer Kane-Mele with real-space sliding hoppings',
      'protocol':'T0 insulating, frozen occupation, static finite-Bz derivative of interband spin-electric susceptibility',
      'exclude':'not nonequilibrium driven DC NMEE; not proven ferroelectric minima; no material claim',
      'physics_relation':'chi_xx(Bz) = -2 sum_occupied_v,empty_c Im[Sx_vc vx_cv]/(Ev-Ec)^2; gamma = finite d/dBz chi',
      'spin_hopping_channels':'full, same-sublattice-only, cross-sublattice-only at IDENTICAL t, SOC, mu, decay; nonadditivity is inclusion-exclusion at fixed parameters',
      'diagnostic_warning':'matrix-element vs denominator split uses instantaneous sorted eigenbasis and is not a universal gauge-invariant separation in degenerate occupied/empty manifolds',
      'base_parameters':asdict(default),
      'coarse_quadrature':meta,'fine_quadrature':meta2,
      'parameter_robustness':results,
      'hopping_channel_competition':channels,
      'weak_hopping_controls':weak,
      'relative_channel_sign_controls':relative_channel_phase,
      'global_layer_gauge_evenness':layer_gauge_evenness,
      'independent_static_kernel_and_quadrature_checks':benchmarks}
    out_path=Path(out_path);out_path.parent.mkdir(parents=True,exist_ok=True)
    out_path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results/sliding_stage4_mechanism.json')
    p.add_argument('--quick',action='store_true');args=p.parse_args();run(args.output,quick=args.quick)
