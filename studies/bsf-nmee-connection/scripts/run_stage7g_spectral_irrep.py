"""Regenerate Stage7g double-valley spectral-sector radial response.

From studies/bsf-nmee-connection:
OPENBLAS_NUM_THREADS=1 python scripts/run_stage7g_spectral_irrep.py
The saved file is recomputed, leaving frozen reference results untouched.
"""
import sys,json,argparse
from pathlib import Path
from dataclasses import replace
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.spectral_irrep_radius import channel_irrep_sectors


def run(quick=False):
    shift=(.4*.23,.4*.41)
    spec=replace(SlidingSpec(),shift=shift,tunneling_eV=.00015625)
    nr,nt=(18,48) if quick else (40,96)
    rows=[]
    for R in (.10,.15,.20,.25,.35,.50,.75,1.):
        x=channel_irrep_sectors(spec,radius=R,nr=nr,nt=nt)
        b=x['channels']['B'];c=x['channels']['C']
        ll=b['per_sector']['low_to_low']['scalar']
        rows.append({'R_inv_nm':R,
                     'B_scalar_full':b['total']['scalar'],
                     'B_scalar_LL':ll,
                     'B_scalar_remote':b['total']['scalar']-ll,
                     'B_four_sectors':{key:val['scalar']
                              for key,val in b['per_sector'].items()},
                     'C_nematic_xx':c['total']['nematic_xx'],
                     'C_four_sectors':{key:val['nematic_xx']
                              for key,val in c['per_sector'].items()},
                     'min_mu_distance_eV':x['min_mu_gap_eV'],
                     'min_group_gap_eV':x['min_remote_group_gap_eV']})
        print(f"R={R:.2f} B={rows[-1]['B_scalar_full']:+.8f} "
              f"LL={ll:+.8f} remote={rows[-1]['B_scalar_remote']:+.8f}",flush=True)
    small=rows[0];large=rows[6]
    ll=large['B_scalar_LL']-small['B_scalar_LL']
    rem=large['B_scalar_remote']-small['B_scalar_remote']
    cancellation=(abs(ll)+abs(rem))/abs(ll+rem)
    return {'protocol':'8-band gapped T=0 equilibrium STATIC spin-electric d_B NMEE',
            'scope':'two K/Kprime disks, not a strict entire-BZ error bound',
            't_eV':spec.tunneling_eV,'u':.4,'nr':nr,'nt':nt,
            'rows':rows,
            'annulus_R0.1_to_R0.75':{'delta_LL':ll,'delta_remote':rem,
                'delta_total':ll+rem,'cancellation_condition_number':cancellation},
            'caveats':'Spectral sectors are isolated-band-group diagnostics; not standalone observables or a dressed four-band theory'}


if __name__=='__main__':
    arg=argparse.ArgumentParser()
    arg.add_argument('--quick',action='store_true')
    args=arg.parse_args()
    result=run(args.quick)
    p=ROOT/'results'/('sliding_stage7g_radial_recomputed_quick.json'
                       if args.quick else 'sliding_stage7g_radial_recomputed.json')
    p.write_text(json.dumps(result,indent=2)+'\n')
    print('saved',p)
