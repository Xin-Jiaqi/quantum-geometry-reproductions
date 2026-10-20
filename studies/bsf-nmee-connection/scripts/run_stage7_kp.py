"""K/K' Taylor models: independent local static-response validation."""
from dataclasses import replace,asdict
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_kp import compare_patch
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.run_stage3_sign_reversal import spec_at


def run(quick=False):
    all_results=[]
    for u in (.1,.4):
      for t in (.005,.02):
        spec=replace(spec_at(u),tunneling_eV=t)
        g=SlidingBilayer(spec).model().reciprocal
        centers=np.array([[1/3,2/3],[2/3,1/3]])@g
        for radius in (.10,.25,.50):
            ks,w,meta=valley_quadrature(g,exterior_n=9,radius=radius,radial_n=16 if quick else 28,angular_n=48 if quick else 72)
            n=meta['n_valley_total']//2
            kk=ks[-2*n:];ww=w[-2*n:]
            for order in (1,2):
                res=compare_patch(spec,kk,ww,centers,n,order=order)
                gamma_full=float(res['full'][0,0]);gamma_kp=float(res['kp'][0,0]);
                trace_full=float(np.trace(res['full']));trace_kp=float(np.trace(res['kp']))
                all_results.append(dict(u=u,t_eV=t,radius_inv_nm=radius,order=order,
                    gamma_full=gamma_full,gamma_kp=gamma_kp,gamma_relative_error=abs(gamma_kp-gamma_full)/max(abs(gamma_full),1e-15),
                    trace_full=trace_full,trace_kp=trace_kp,
                    trace_absolute_error=abs(trace_kp-trace_full),
                    min_mu_gap_eV=res['min_mu_gap'].tolist(),n_points=2*n))
                print(f'u {u} t {t} r {radius} order {order} exact {gamma_full:.7e} kp {gamma_kp:.7e} err {all_results[-1]["gamma_relative_error"]:.3g}',flush=True)
    return all_results

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');parser.add_argument('--output',type=Path,default=ROOT/'results/sliding_stage7_valley_kp.json');args=parser.parse_args()
    blob={'scope':'zero-K gapped static response, local K/K prime disk, no driven-DC claim','data':run(quick=args.quick)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(blob,indent=2)+'\n')
