"""Re-run the Stage7d valley-vertex audit and hopping-ratio sign test.

From studies/bsf-nmee-connection: python scripts/run_stage7d_vertex_ratio.py
This is two-valley DISK static NMEE, not full-BZ or driven DC.
"""
from pathlib import Path
from dataclasses import replace
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingSpec, SlidingBilayer
from bsf_nmee.valley_vertex_symmetry import (
    chiral_valley_summary,rashba_spin_rotation_check)
from bsf_nmee.channel_spectral_mechanism import disk_channel_sectors


def run():
    aa=replace(SlidingSpec(),shift=(0.,0.),tunneling_eV=.005)
    rows=chiral_valley_summary(aa)
    shifted=replace(SlidingSpec(),shift=(.4*.23,.4*.41),
                    tunneling_eV=.0025)
    k=np.array([1/3,2/3])@SlidingBilayer(shifted).model().reciprocal+np.array([.13,.07])
    spin_check=rashba_spin_rotation_check(shifted,k)
    scans=[]
    for ratio in (.50,.53,.54,.57):
        r=disk_channel_sectors(shifted,cross_scale=ratio,radius=.75,
                              radial_n=26,angular_n=72)
        scans.append(dict(ratio=ratio,gamma=r['total'],
                          gamma_over_s2=r['total']/.0025**2,
                          mu_gap_eV=r['sampled_min_mu_distance_eV'],
                          sectors=r['sectors']))
    return {'protocol':'T=0 gapped equilibrium static two K/Kprime disks',
            'valley_reference':rows,'rashba_unitary':spin_check,
            'u':.4,'s_eV':.0025,'ratio_scan':scans,
            'caution':'This is an artificial orbital-channel amplitude comparison, not a real device or a driven-DC coefficient.'}


if __name__=='__main__':
    out=run()
    target=ROOT/'results/sliding_stage7d_vertex_ratio.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('saved',target)
    for r in out['ratio_scan']:
        print('c/s',r['ratio'],'gamma',r['gamma'],'gap',r['mu_gap_eV'])
