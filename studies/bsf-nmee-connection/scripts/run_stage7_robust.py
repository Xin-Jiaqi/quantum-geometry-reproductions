"""Small symmetry-preserving model perturbation screen (full Brillouin zone)."""
from dataclasses import replace
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.symmetry_scaling import static_inplane_matrix
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.run_stage3_sign_reversal import spec_at

def run():
    base=spec_at(.1);g=SlidingBilayer(base).model().reciprocal
    k,w,meta=valley_quadrature(g,exterior_n=41,radius=1.5,radial_n=35,angular_n=64)
    knobs={'reference':{},'rashba_minus20':{'rashba_soc_eV':.016},
           'rashba_plus20':{'rashba_soc_eV':.024},
           'intrinsic_minus20':{'intrinsic_soc_eV':.008},
           'intrinsic_plus20':{'intrinsic_soc_eV':.012},
           'decay_minus12p5':{'decay_nm':.14},
           'decay_plus12p5':{'decay_nm':.18}}
    results={}
    for name,config in knobs.items():
        dat=[]
        for u in (.1,.4):
            sp=replace(spec_at(u),tunneling_eV=.005,**config)
            try:
                val=static_inplane_matrix(sp,k,w)
                dat.append({'u':u,'gamma_xz_x':float(val['gamma'][0,0]),'trace':float(np.trace(val['gamma'])),
                            'min_mu_distance_eV':val['mu_distance']})
            except Exception as e:
                dat.append({'u':u,'error':repr(e)})
        results[name]={'parameter_updates':config,'data':dat}
        print(name,[(v['u'],v.get('gamma_xz_x'),v.get('error')) for v in dat],flush=True)
    return {'protocol':'zero-temperature insulating full-BZ static finite-B NMEE','integration':meta,'values':results}
if __name__=='__main__':
 x=run();(ROOT/'results/sliding_stage7_robust.json').write_text(json.dumps(x,indent=2)+'\n')
