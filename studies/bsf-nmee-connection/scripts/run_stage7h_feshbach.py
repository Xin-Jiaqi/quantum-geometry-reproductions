"""Regenerate the Stage7h local exact Feshbach source-response checks.

Run inside studies/bsf-nmee-connection:
OPENBLAS_NUM_THREADS=1 python scripts/run_stage7h_feshbach.py
Does not overwrite the frozen Stage7h benchmark JSON.
"""
from pathlib import Path
from dataclasses import replace
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from bsf_nmee.feshbach_response import feshbach_static_at_k


def run():
    rows=[]
    for t in (.005,.020):
        spec=replace(SlidingSpec(),shift=(.4*.23,.4*.41),
                     tunneling_eV=t)
        g=SlidingBilayer(spec).model().reciprocal
        centers=np.array([[1/3,2/3],[2/3,1/3]])@g
        for nu,center in enumerate(centers):
            for q in ((0.,0.),(.05,.02),(.25,0.)):
                val=feshbach_static_at_k(
                    spec,center+q,step_T=.25,contour_n=64)
                val['valley']=nu;val['q_inv_nm']=list(q)
                rows.append(val)
                print('t',t,'nu',nu,'q',q,'full',val['gamma_eigen'],
                      'dressed-only',val['gamma_dressed_only_diagnostic'],
                      'contacts',val['gamma_missing_contacts'],flush=True)
    return {'scope':'gapped eight-band equilibrium STATIC finite-B NMEE, local K/Kprime points',
            'source':'exact Schur identity and independently diagonalized 8band Kubo',
            'rows':rows,'caveats':['not a complete four-band response theory',
                'dressed-only is an algebraic diagnostic, not observable',
                'not integrated full BZ and not driven DC']}


if __name__=='__main__':
    data=run()
    p=ROOT/'results/sliding_stage7h_feshbach_recomputed.json'
    p.write_text(json.dumps(data,indent=2)+'\n')
    print('saved',p)
