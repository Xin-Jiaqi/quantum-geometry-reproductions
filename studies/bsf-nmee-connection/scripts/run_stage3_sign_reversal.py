"""Stage 3: insulating sliding-path reversal and orbital-channel controls.

Only a *static insulating* finite-B derivative of the spin-electric interband
susceptibility is considered. Not the unsettled driven DC PRL NMEE coefficient.
All shifts are fractional and all observed quantities retain stage-2 units.
"""
from __future__ import annotations
import argparse
import json
from dataclasses import replace
from datetime import date
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer, SlidingSpec
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.sliding_static_scan import response_weighted

SOURCE_SHA = 'c429336f1a3e4827559c89a75425ced7af69c696'


def spec_at(u: float, *, same_only=False, tunneling=0.020) -> SlidingSpec:
    return SlidingSpec(shift=(0.23*u, 0.41*u), tunneling_eV=tunneling,
                       cross_sublattice=not same_only)


def hopping_norm_ratio(spec: SlidingSpec, n=41) -> dict:
    """Match BZ average Tr(T^dagger T), not a single-point hopping amplitude."""
    obj = SlidingBilayer(spec)
    frac = (np.arange(n)+0.5)/n
    k = np.stack(np.meshgrid(frac,frac,indexing='ij'),axis=-1).reshape(-1,2) @ obj.model().reciprocal
    t = obj.tunneling(k)
    same = t.copy()
    same[:,:2,2:] = 0
    same[:,2:,:2] = 0
    p_full=float(np.mean(np.sum(np.abs(t)**2,axis=(1,2))))
    p_same=float(np.mean(np.sum(np.abs(same)**2,axis=(1,2))))
    if not (p_full>p_same>0):
        raise AssertionError('Inconsistent interlayer-hopping norm decomposition')
    return dict(bz_mesh=n, total_norm_sq_eV2=p_full, same_norm_sq_eV2=p_same,
                cross_norm_sq_eV2=p_full-p_same,
                same_amplitude_rescale=float(np.sqrt(p_full/p_same)))


def integrate_partitioned(spec,k,w,meta,step=.25):
    ext=meta['n_exterior_kept']; valley=meta['n_valley_total']//2
    parts={}; j=0
    for name,n in [('exterior',ext),('K',valley),('Kprime',valley)]:
        r=response_weighted(spec,k[j:j+n],w[j:j+n],step=step)
        parts[name]=float(r['response']['total_x']); j+=n
    if j!=len(k):raise AssertionError('quadrature partition mismatch')
    return parts


def run(output:Path,quick=False):
    base=SlidingBilayer().model()
    if quick:
        coarse=(41,35,64,1.5)
        fine=(61,60,80,1.5)
        steps=[0,.1,.2,.4,1.]
        root_steps=4
    else:
        coarse=(41,35,64,1.5)
        fine=(61,60,80,1.5)
        steps=[0,.05,.1,.15,.175,.2,.225,.25,.3,.35,.4,.5,.6,.75,1.,1.25]
        root_steps=11
    k,w,meta=valley_quadrature(base.reciprocal,coarse[0],coarse[3],coarse[1],coarse[2])
    k2,w2,meta2=valley_quadrature(base.reciprocal,fine[0],fine[3],fine[1],fine[2])
    for x in (meta,meta2):
        if abs(x['normalization_error'])/x['bz_area_weight']>3e-3:
            raise AssertionError('BZ quadrature normalization inaccurate')
    cache={}
    def eval_u(u,points=k,weights=w,step=.25):
        key=(round(float(u),12),len(points),step)
        if key not in cache:
            cache[key]=response_weighted(spec_at(u),points,weights,step=step)
        return cache[key]
    scan=[]
    for u in steps:
        result=eval_u(u)
        if result['min_mu_to_band_eV']<=0:
            raise AssertionError(f'Chemical potential outside gap for u={u}')
        print('slide',u,format(result['response']['total_x'],'.9e'),flush=True)
        scan.append({'path_parameter':u,'shift_fractional':list(spec_at(u).shift),
                     'gamma':float(result['response']['total_x']),
                     'top':float(result['response']['top_x']),
                     'bottom':float(result['response']['bottom_x']),
                     'sampled_direct_gap_eV':float(result['sampled_middle_band_gap_eV']),
                     'sampled_min_mu_to_band_eV':float(result['min_mu_to_band_eV'])})
    left=.15;right=.20
    fl=eval_u(left)['response']['total_x'];fr=eval_u(right)['response']['total_x']
    if fl>=0 or fr<=0:raise AssertionError('no sign change in expected root bracket')
    for _ in range(root_steps):
        mid=(left+right)/2
        fm=eval_u(mid)['response']['total_x']
        print('bisect',mid,format(fm,'.9e'),flush=True)
        if fm>0:right=mid
        else:left=mid
    center=(left+right)/2
    root_check=eval_u(center)
    sign_nodes={'root_bracket':[left,right],'estimated_zero_u':center,
                'residual_gamma_at_center':float(root_check['response']['total_x']),
                'sampled_mu_gap_eV':float(root_check['min_mu_to_band_eV'])}
    convergence={};partitions={};step_control={}
    for u in (.1,.2,.4,1.):
        low=eval_u(u);high=eval_u(u,k2,w2)
        a=low['response']['total_x'];b=high['response']['total_x']
        print('convergence',u,format(a,'.9e'),format(b,'.9e'),flush=True)
        convergence[str(u)]={'coarse_gamma':float(a),'fine_gamma':float(b),
                             'absolute_change':float(abs(a-b)),
                             'relative_change_if_nonzero':float(abs(a-b)/max(abs(b),1e-12)),
                             'fine_sampled_mu_gap_eV':float(high['min_mu_to_band_eV'])}
        partitions[str(u)]=integrate_partitioned(spec_at(u),k,w,meta)
    for u in (.1,.2,.4):
        a=eval_u(u)['response']['total_x'];b=eval_u(u,k,w,step=.125)['response']['total_x']
        step_control[str(u)]={'step_0p25_T':float(a),'step_0p125_T':float(b),
                               'absolute_change':float(abs(a-b))}
    orb={}
    for u in (.4,1.):
        sp=spec_at(u)
        ratio=hopping_norm_ratio(sp)
        full=eval_u(u)['response']['total_x']
        same_raw=response_weighted(spec_at(u,same_only=True),k,w)
        same_matched=response_weighted(spec_at(u,same_only=True,tunneling=sp.tunneling_eV*ratio['same_amplitude_rescale']),k,w)
        if same_matched['min_mu_to_band_eV']<=0:
            raise AssertionError('Norm-matched same-sublattice control is metallic')
        orb[str(u)]={'normalization':ratio, 'full_gamma':float(full),
                     'same_only_raw_gamma':float(same_raw['response']['total_x']),
                     'same_only_norm_matched_gamma':float(same_matched['response']['total_x']),
                     'matched_same_tunneling_eV':float(sp.tunneling_eV*ratio['same_amplitude_rescale']),
                     'matched_same_sampled_mu_gap_eV':float(same_matched['min_mu_to_band_eV'])}
    parity={}
    for u in (.1,.4,1.):
        positive=eval_u(u)
        neg=response_weighted(spec_at(-u),k,w)
        parity[str(u)]={'positive':float(positive['response']['total_x']),
                        'negative':float(neg['response']['total_x']),
                        'odd_residual':float(positive['response']['total_x']+neg['response']['total_x']),
                        'top_partner_error':float(positive['response']['top_x']+neg['response']['bottom_x'])}
    report={'protocol':'T=0 insulating interband E-linear susceptibility, static Bz finite difference, frozen valence occupancy',
            'NOT':'author PRL finite-temperature driven nonlinear DC response',
            'source_parent_github_sha':SOURCE_SHA,'date':str(date.today()),
            'fixed_controls':{'tunneling_eV':.020,'layer_potential_difference_eV':0.,
                              'mu_eV':spec_at(0).mu_eV,'B_step_T':.25,
                              'trajectory':'tau(u)=(0.23u,0.41u) in fractional primitive lattice',
                              'unit':'(nm V T)^-1 under current NMEE Pauli-spin conventions'},
            'quadratures':{'coarse':dict(zip(['exterior_n','radial_n','angular_n','radius_inv_nm'],coarse))|meta,
                            'fine':dict(zip(['exterior_n','radial_n','angular_n','radius_inv_nm'],fine))|meta2},
            'scan':scan,'zero_crossing':sign_nodes,'mesh_convergence':convergence,
            'magnetic_field_difference':step_control,'valley_partition':partitions,
            'equal_hopping_norm_controls':orb,'mirror_oddness':parity}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(output),'zero_u':center,'mesh_convergence':convergence,
                      'orbital_controls':orb,'mirror':parity},indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--quick',action='store_true',help='faster but not publication-quality smoke benchmark')
    p.add_argument('--output',type=Path,default=ROOT/'results/sliding_stage3_sign_reversal.json')
    options=p.parse_args()
    run(options.output,quick=options.quick)
