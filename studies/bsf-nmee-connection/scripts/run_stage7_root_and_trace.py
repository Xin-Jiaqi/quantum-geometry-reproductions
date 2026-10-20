"""Local k.p predicts sliding-root location and C3 trace/covariance scaling."""
from dataclasses import replace
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.valley_kp import compare_patch
from scripts.run_stage3_sign_reversal import spec_at


def run():
    base=SlidingBilayer();g=base.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@g
    k,w,m=valley_quadrature(g,exterior_n=7,radius=.75,radial_n=28,angular_n=72)
    n=m['n_valley_total']//2; k=k[-2*n:];w=w[-2*n:]
    root_data=[]
    for u in (.15,.175,.2,.225):
        s=spec_at(u)
        r=compare_patch(s,k,w,centers,n,order=2)
        root_data.append({'u':u,'full_disks_xx':float(r['full'][0,0]),'kp2_disks_xx':float(r['kp'][0,0])})
        print('root sample',u,root_data[-1],flush=True)
    a=.175;b=.2
    for _ in range(8):
        mid=(a+b)/2
        s=spec_at(mid);r=compare_patch(s,k,w,centers,n,order=2)
        if r['kp'][0,0]>0:b=mid
        else:a=mid
    root=.5*(a+b)
    trace=[]
    for u in (.025,.05,.1):
        s=replace(spec_at(u),tunneling_eV=.005)
        r=compare_patch(s,k,w,centers,n,order=2)
        def nem(m):
            sym=(m+m.T)/2; tr=np.trace(m)/2
            return float(np.linalg.norm(sym-tr*np.eye(2)))
        trace.append({'u':u,'full_disk_trace':float(np.trace(r['full'])),'kp2_disk_trace':float(np.trace(r['kp'])),
                      'full_disk_nematic':nem(r['full']),'kp2_disk_nematic':nem(r['kp'])})
        print('trace',u,trace[-1],flush=True)
    return {'path_zero_kp2':[a,b],'estimated_zero_kp2':root,'original_stage3_BZ_root_approx':.1974,
            'root_scan':root_data,'trace_scan':trace}
if __name__=='__main__':
 x=run();(ROOT/'results/sliding_stage7_root_and_trace.json').write_text(json.dumps(x,indent=2)+'\n')
