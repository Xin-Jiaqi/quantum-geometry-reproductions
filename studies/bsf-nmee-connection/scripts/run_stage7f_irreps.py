"""Regenerate Stage7f channel tensors and C3/mirror physics checks.

Run from studies/bsf-nmee-connection:
 OPENBLAS_NUM_THREADS=1 python scripts/run_stage7f_irreps.py --quick
 OPENBLAS_NUM_THREADS=1 python scripts/run_stage7f_irreps.py
Will not overwrite the frozen reference results.
"""
from dataclasses import replace
from pathlib import Path
import argparse,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.irrep_channel_tensor import (
    analytic_tensor_ABC,rotate_slide,tensor_covariance_residual)


def spec_for(u):
    return replace(SlidingSpec(),shift=(.23*u,.41*u))


def run(quick=False):
    n,t=(12,36) if quick else (26,72)
    series=[]
    for u in (.02,.05,.1,.2,.4):
        response=analytic_tensor_ABC(spec_for(u),nr=n,nt=t)
        tensor=response["tensors"]
        row={"u":u,"data":tensor,
             "coefficient_Gamma_xz_x_at_same_equals_cross":sum(
                  tensor[x]["matrix"][0][0] for x in ("A","B","C"))}
        series.append(row)
        print("u",u,"Gamma_xz;x /t²",row["coefficient_Gamma_xz_x_at_same_equals_cross"],flush=True)
    reference=analytic_tensor_ABC(spec_for(.1),nr=n,nt=t)
    rotated=analytic_tensor_ABC(rotate_slide(spec_for(.1),120),nr=n,nt=t)
    rotated60=analytic_tensor_ABC(rotate_slide(spec_for(.1),60),nr=n,nt=t)
    reversed_=analytic_tensor_ABC(spec_for(-.1),nr=n,nt=t)
    tests={
        "C3_120_max_error":tensor_covariance_residual(reference,rotated,120),
        "mirror_reverse_max_error":tensor_covariance_residual(reference,reversed_,0,odd=True),
        "60_equivalent_240_and_mirror_max_error":tensor_covariance_residual(reference,rotated60,240,odd=True),
    }
    print("symmetry checks",tests)
    return {"protocol":"synthetic 8-band, T=0 gapped equilibrium static derivative, two K/Kprime disks",
            "radius_inv_nm":.75,"nr":n,"nt":t,"series":series,"checks":tests,
            "caveat":"channel basis is fixed A/B; not driven-DC and not verified ferroelectric energetics"}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('--quick',action='store_true')
    args=parser.parse_args()
    result=run(args.quick)
    name='sliding_stage7f_recomputed_quick.json' if args.quick else 'sliding_stage7f_recomputed.json'
    dest=ROOT/'results'/name
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print("saved",dest)
