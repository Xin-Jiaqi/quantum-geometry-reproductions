"""Weak-coupling mixed-channel derivatives of the *static* gapped response.

C_sc = lim_{t->0} [Gamma(t_s=t,t_c=t)-Gamma(t_s=t,t_c=-t)]/(2*t**2)
by the exact layer-sign unitary Gamma(-s,-c)=Gamma(s,c).  The odd/even
channel parity is a controlled Hamiltonian comparison, not a basis-independent
assignment to isolated bands.  No assertion regarding the PRL driven DC NMEE.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import replace, asdict
import argparse, json, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.channel_perturbation import evaluate_scaled_channels
from scripts.run_stage3_sign_reversal import spec_at


def run(out, *, quick=False):
    g=SlidingBilayer().model().reciprocal
    coarse=(29,1.5,25,48) if quick else (41,1.5,35,64)
    fine=(41,1.5,35,64) if quick else (61,1.5,60,80)
    k,w,meta=valley_quadrature(g,*coarse)
    kf,wf,fmeta=valley_quadrature(g,*fine)
    t_values=[.0025,.005,.01,.02] if quick else [.00125,.0025,.005,.01,.02]
    us=[.1,.2,.4]
    data={}
    for u in us:
        rows=[]
        for t in t_values:
            spec=replace(spec_at(u),tunneling_eV=t)
            pos=evaluate_scaled_channels(spec,k,w)
            neg=evaluate_scaled_channels(spec,k,w,cross_scale=-1)
            if quick and t not in (.0025,.02):
                so=co=None
            else:
                so=evaluate_scaled_channels(spec,k,w,cross_scale=0)['gamma']
                co=evaluate_scaled_channels(spec,k,w,same_scale=0)['gamma']
            gp,gm=pos['gamma'],neg['gamma']
            even=.5*(gp+gm);odd=.5*(gp-gm)
            rows.append({'t_eV':t,'plus':gp,'cross_sign_reversed':gm,
                         'even_under_cross_sign':even, 'odd_under_cross_sign':odd,
                         'total_coefficient_over_t_sq':gp/t**2,
                         'even_coefficient_over_t_sq':even/t**2,
                         'mixed_coefficient_over_t_sq':odd/t**2,
                         'same_only':so,'cross_only':co,
                         'nonadditivity':None if so is None else gp-so-co,
                         'min_mu_gap_eV':min(pos['sampled_min_mu_distance_eV'],neg['sampled_min_mu_distance_eV'])})
        # Fit leading Taylor coefficients C_2 with linear dependence on t^2;
        # numerical regression is purely an extrapolation, not proof of analyticity.
        fitted={}
        for label in ('total_coefficient_over_t_sq','even_coefficient_over_t_sq','mixed_coefficient_over_t_sq'):
            m=3 if len(rows)>3 else 2
            xx=np.array([r['t_eV']**2 for r in rows[:m]])
            yy=np.array([r[label] for r in rows[:m]])
            slope,intercept=np.polyfit(xx,yy,deg=1)
            fitted[label.replace('_over_t_sq','')+'_t_to_zero']={'intercept':float(intercept),'quartic_slope':float(slope),
                        'raw_first':float(yy[0]),'fit_points_eV':[r['t_eV'] for r in rows[:m]]}
        data[f'{u:.2f}']={'sweep':rows,'taylor_extrapolation':fitted}
        print('u',u,'t0 total',fitted['total_coefficient_t_to_zero']['intercept'],
              't0 odd-mixed',fitted['mixed_coefficient_t_to_zero']['intercept'],flush=True)
    # Gauge parity + numerical quadrature + Zeeman difference convergence checks
    checks=[]
    for u,t in ((.1,.005),(.4,.005),(.1,.02),(.4,.02)):
        s=replace(spec_at(u),tunneling_eV=t)
        # smaller B field step probes B slope, finer quadrature probes BZ integration.
        p=evaluate_scaled_channels(s,k,w)['gamma']
        pf=evaluate_scaled_channels(s,kf,wf)['gamma']
        bhalf=evaluate_scaled_channels(s,k,w,step_T=.125)['gamma']
        checks.append({'u':u,'t_eV':t,'coarse':p,'fine':pf,'fine_minus_coarse':pf-p,
                       'half_B_step':bhalf,'half_B_step_minus_default':bhalf-p})
        print('check',u,t,'mesh delta',pf-p,'B delta',bhalf-p,flush=True)
    # Four quadrants confirm mixed finite-difference identity by independent numbers.
    s=replace(spec_at(.1),tunneling_eV=.005)
    vals={f'{i:+d}{j:+d}':evaluate_scaled_channels(s,k,w,same_scale=i,
                                                cross_scale=j)['gamma']
          for i,j in [(1,1),(1,-1),(-1,1),(-1,-1)]}
    sc4=(vals['+1+1']-vals['+1-1']-vals['-1+1']+vals['-1-1'])/(4*s.tunneling_eV**2)
    assert abs(vals['+1+1']-vals['-1-1'])<1e-11 and abs(vals['+1-1']-vals['-1+1'])<1e-11
    assert abs(sc4-(vals['+1+1']-vals['+1-1'])/(2*s.tunneling_eV**2))<1e-8
    # Expected phenomenology at finite coupling: odd is negative at BOTH sliding endpoints;
    # the even contribution overtakes it by u=0.4.
    p1=data['0.10']['sweep'][-1];p4=data['0.40']['sweep'][-1]
    assert p1['odd_under_cross_sign']<0 and p4['odd_under_cross_sign']<0
    assert p1['even_under_cross_sign']>0 and p4['even_under_cross_sign']>0
    assert p1['plus']<0<p4['plus']
    out=Path(out);out.parent.mkdir(exist_ok=True,parents=True)
    blob={'scope':'toy 8band Kane-Mele; 0K insulator; equilibrium/static mixed curvature; finite field difference',
          'excluded':'not author driven DC NMEE / not proven ferroelectric / not a material prediction',
          'conventions':'s and c multiply same- and cross-sublattice interlayer amplitudes, respectively; basis fixed; u is fraction of (0.23,0.41)',
          'base_spec':asdict(spec_at(.1)),'coarse':{'params':coarse,'meta':meta},'fine':{'params':fine,'meta':fmeta},
          'data':data,'convergence':checks,'four_quadrant':{'responses':vals,'finite_mixed_coefficient':sc4}}
    out.write_text(json.dumps(blob,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return blob

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results/sliding_stage5_weak_coupling.json');p.add_argument('--quick',action='store_true');args=p.parse_args()
    run(args.output,quick=args.quick)
