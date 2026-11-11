"""Stage 7i: reproduce adaptive K/K' Wilson/Kubo static tensor checks.

From studies/bsf-nmee-connection:
 OPENBLAS_NUM_THREADS=1 python scripts/run_stage7i_wilson.py
Writes only results/sliding_stage7i_recomputed.json.
"""
from pathlib import Path
from dataclasses import replace
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.mixed_wilson_geometry import compare_valley_wilson_kubo


def run():
    results=[]
    for u,t in ((.1,.005),(.4,.005),(.4,.02)):
        spec=replace(SlidingSpec(),shift=(u*.23,u*.41),tunneling_eV=t)
        for dk,dl in ((.01,.00025),(.005,.000125),(.0025,.0000625)):
            out=compare_valley_wilson_kubo(
                spec,dk=dk,dl=dl,radius=.5,radial_n=22,angular_n=60)
            results.append(out)
            print('u',u,'t',t,'dk',dk,'Wilson',out['wilson'],
                  'Kubo',out['kubo'],'error',out['relative_difference'],flush=True)
    matrix_checks=[]
    for u in (.1,.4):
        spec=replace(SlidingSpec(),shift=(u*.23,u*.41),tunneling_eV=.005)
        rows=[]
        for spin in ('x','y'):
            rows.append([compare_valley_wilson_kubo(
                spec,spin_axis=spin,k_axis=j,
                dk=.0025,dl=.0000625,radius=.5,radial_n=22,angular_n=60)
                for j in (0,1)])
        matrix_checks.append({'u':u,'tensor':rows})
    return {'protocol':'gapped T=0 equilibrium static mixed Berry curvature',
            'domain':'K and Kprime radius 0.5 nm^-1, NOT entire Brillouin zone',
            'refinements':results,'tensor_checks':matrix_checks,
            'DFT_status':'No actual material or first-principles calculation has been run'}


if __name__=='__main__':
    data=run()
    target=Path(__file__).resolve().parents[1]/'results/sliding_stage7i_recomputed.json'
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('saved',target)
