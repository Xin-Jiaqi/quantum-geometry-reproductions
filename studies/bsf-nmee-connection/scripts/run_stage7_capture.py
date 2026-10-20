"""Valley-capture fraction: gapped static NMEE (no published driven-DC claim)."""
from dataclasses import replace
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.symmetry_scaling import static_inplane_matrix
from bsf_nmee.valley_kp import compare_patch
from scripts.run_stage3_sign_reversal import spec_at


def run():
    out=[]
    for u in (.1,.4):
        for t in (.005,.02):
            sp=replace(spec_at(u),tunneling_eV=t)
            obj=SlidingBilayer(sp); g=obj.model().reciprocal
            centers=np.array([[1/3,2/3],[2/3,1/3]])@g
            all_k,w,meta=valley_quadrature(g,exterior_n=61,radius=1.5,radial_n=60,angular_n=80)
            full=static_inplane_matrix(sp,all_k,w)['gamma']
            for radius in (.1,.25,.5,.75):
                kk,ww,m=valley_quadrature(g,exterior_n=9,radius=radius,radial_n=28,angular_n=72)
                n=m['n_valley_total']//2
                r=compare_patch(sp,kk[-2*n:],ww[-2*n:],centers,n,order=2)
                d={"u":u,"t_eV":t,"radius_inv_nm":radius,
                   "BZ_gamma_xx":float(full[0,0]),"BZ_trace":float(np.trace(full)),
                   "disk_exact_gamma_xx":float(r['full'][0,0]),"disk_kp2_gamma_xx":float(r['kp'][0,0]),
                   "disk_exact_trace":float(np.trace(r['full'])),"disk_kp2_trace":float(np.trace(r['kp'])),
                   "disk_fraction_signed_xx":float(r['full'][0,0]/full[0,0]),
                   "disk_fraction_signed_trace":float(np.trace(r['full'])/np.trace(full)) if abs(np.trace(full))>1e-15 else None,
                   "kp2_vs_disk_relative_error_xx":float(abs(r['kp'][0,0]-r['full'][0,0])/abs(r['full'][0,0])),
                   "kp2_vs_BZ_relative_error_xx":float(abs(r['kp'][0,0]-full[0,0])/abs(full[0,0])),
                   "min_mu_distance_eV":r['min_mu_gap'].tolist()}
                out.append(d)
                print('u',u,'t',t,'r',radius,'BZ',f'{full[0,0]:.6g}','disk',f'{r["full"][0,0]:.6g}','kp2',f'{r["kp"][0,0]:.6g}','ratio',f'{d["disk_fraction_signed_xx"]:.4f}',flush=True)
    return out
if __name__=='__main__':
    out=run();p=ROOT/'results/sliding_stage7_capture.json';p.write_text(json.dumps(out,indent=2)+'\n');print('saved',p)
