"""Independent scalar benchmark and scope-controlled symmetry/hybridization audit."""
from __future__ import annotations
import sys,json,time,hashlib,datetime
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from bsf_nmee.sliding_tight_binding import SlidingSpec, SlidingBilayer, M_GLOBAL
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.sliding_static_scan import response_weighted


def main():
    G=SlidingBilayer().model().reciprocal
    configs=[('base',.028,True,(.23,.41),.25),('decoupled',0.,True,(.23,.41),.25),
             ('same_sublattice',.028,False,(.23,.41),.25),
             ('weak10',.010,True,(.23,.41),.25),('weak20',.020,True,(.23,.41),.25),
             ('strong40',.040,True,(.23,.41),.25),
             ('mirror_partner',.028,True,(-.23,-.41),.25),
             ('step_0p125',.028,True,(.23,.41),.125),
             ('small05',.005,True,(.23,.41),.25)]
    meshes=[('coarse',81,90,120,1.5)]
    results={'models':{},'convergence':{},'quadrature':{},'physics_scope':'T=0 insulating equilibrium static finite-B derivative (frozen occupied projector), NOT driven NMEE DC'}
    start=time.time()
    for label,n,nr,na,R in meshes:
        pts,w,meta=valley_quadrature(G,n,R,nr,na)
        result=response_weighted(SlidingSpec(),pts,w)
        results['convergence'][label]=result
        results['quadrature'][label]={**meta,'n_exterior':n,'n_radial':nr,'n_angular':na,'radius_1_over_nm':R}
        print('mesh',label,result['response']['total_x'],flush=True)
        if label!='coarse':continue
        for name,t,cross,shift,bstep in configs:
            spec=SlidingSpec(tunneling_eV=t,cross_sublattice=cross,shift=shift)
            result=response_weighted(spec,pts,w,step=bstep)
            results['models'][name]={**result,'inputs':{'tunneling_eV':t,'cross_sublattice':cross,'shift':shift,'B_step_T':bstep}}
            print('model',name,result['response']['total_x'],flush=True)
    s=SlidingBilayer(SlidingSpec()); m=SlidingBilayer(SlidingSpec(shift=(-.23,-.41)))
    rng=np.random.default_rng(12345)
    k=rng.uniform(-10,10,size=(10,2))
    mirror_diff=np.max(abs(m.matrices(k)-M_GLOBAL@s.matrices(k)@M_GLOBAL.conj().T))
    mirror_vel=np.max(abs(m.matrices(k,0)-M_GLOBAL@s.matrices(k,0)@M_GLOBAL.conj().T))
    shifted=SlidingBilayer(SlidingSpec(shift=(1.23,.41)))
    periodic_diff=np.max(abs(shifted.matrices(k)-s.matrices(k)))
    periodic_vdiff=np.max(abs(shifted.matrices(k,0)-s.matrices(k,0)))
    dk=np.zeros_like(k);dk[:,0]=1e-5
    speed_diff=np.max(abs((s.matrices(k+dk)-s.matrices(k-dk))/(2e-5)-s.matrices(k,0)))
    results['exact_checks']={'mirror_hamiltonian_eV':float(mirror_diff),'mirror_velocity_eV_nm':float(mirror_vel),
        'sliding_plus_a1_H_eV':float(periodic_diff),'sliding_plus_a1_velocity_eV_nm':float(periodic_vdiff),
        'velocity_central_difference_error_eV_nm':float(speed_diff)}
    results['elapsed_s']=round(time.time()-start,1)
    path=ROOT/'results/sliding_stage2_validated.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(results,indent=2)+'\n')
    return results
if __name__=='__main__':main()
