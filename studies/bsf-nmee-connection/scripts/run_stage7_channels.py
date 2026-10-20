"""Nonadditive AA/BB--AB/BA static coefficient in K/K' patches vs full BZ."""
from dataclasses import replace
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_kp import compare_patch
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.channel_perturbation import evaluate_scaled_channels
from scripts.run_stage3_sign_reversal import spec_at

def run():
    rows=[]
    g=SlidingBilayer().model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@g
    kk,ww,meta=valley_quadrature(g,exterior_n=9,radius=.75,radial_n=28,angular_n=72)
    n=meta['n_valley_total']//2; kk=kk[-2*n:];ww=ww[-2*n:]
    K,W,_=valley_quadrature(g,exterior_n=41,radius=1.5,radial_n=35,angular_n=64)
    for u in (.1,.2,.4):
        s=replace(spec_at(u),tunneling_eV=.005)
        cases={}
        for name,sc,cc in [('full',1,1),('cross_reverse',1,-1),('same_only',1,0),('cross_only',0,1),('zero',0,0)]:
            r=compare_patch(s,kk,ww,centers,n,order=2,same_scale=sc,cross_scale=cc)
            full_bz=evaluate_scaled_channels(s,K,W,same_scale=sc,cross_scale=cc)['gamma']
            cases[name]={'BZ':full_bz,'disk_full':float(r['full'][0,0]),'disk_kp2':float(r['kp'][0,0]),'mu_gap_eV':r['min_mu_gap'].tolist()}
        def odd(which):return (cases['full'][which]-cases['cross_reverse'][which])/2
        def even(which):return (cases['full'][which]+cases['cross_reverse'][which])/2
        entry={'u':u,'t_eV':s.tunneling_eV,'radius_inv_nm':.75,'cases':cases}
        for kind in ('BZ','disk_full','disk_kp2'):
            entry[kind+'_mixed_coeff']=float(odd(kind)/s.tunneling_eV**2)
            entry[kind+'_even_coeff']=float(even(kind)/s.tunneling_eV**2)
            entry[kind+'_nonadditivity']=float(cases['full'][kind]-cases['same_only'][kind]-cases['cross_only'][kind]+cases['zero'][kind])
        rows.append(entry)
        print('u',u,'BZ mixed C',f'{entry["BZ_mixed_coeff"]:.6g}','kp2 mixed C',f'{entry["disk_kp2_mixed_coeff"]:.6g}',
              'BZ even C',f'{entry["BZ_even_coeff"]:.6g}','kp2 even C',f'{entry["disk_kp2_even_coeff"]:.6g}',flush=True)
    return rows
if __name__=='__main__':
    x=run();(ROOT/'results/sliding_stage7_channels.json').write_text(json.dumps(x,indent=2)+'\n')
