"""Independent PT-6-band finite-eta Laurent JSON benchmark."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from methods.pt_synthetic_v31 import PTThreeDoublets
from methods.laurent_v31 import (adiabatic_laurent,explicit_three_cluster_term,
    independent_population_delta,second_order_outer_rate,reconstruct_response)

def run():
    f=PTThreeDoublets()
    rng=np.random.default_rng(14)
    a=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    h=(a+a.conj().T)/2;O=(h+f.J@h.conj()@f.J.conj().T)/2
    r=adiabatic_laurent(f,.37,0.,mu=.6,temperature=.32)
    trace=lambda q:float(np.trace(O@q).real)
    data={'paper':'Xiang/Jin/Wang PRL 135 146303 (2025)',
        'model':'synthetic exact PT six-band three doublets, NOT MnBi2Te4',
        'protocol':r['protocol'],'scientific_scope':r['scientific_scope'],
        'local_k':.37,'mu':.6,'temperature':.32,'cluster_ranks':r['cluster_ranks'],
        'local':{q:trace(r[q]) for q in ('A','B','C_adiabatic','C_half_population_candidate','delta_printed_candidate')},
        'three_cluster_norm':float(np.max(abs(explicit_three_cluster_term(f,.37,0.,mu=.6,temperature=.32)))),
        'population_independent_residual':float(np.max(abs(r['delta_printed_candidate']-
            independent_population_delta(f,.37,0.,mu=.6,temperature=.32))))}
    data['finite_eta_to_constant']=[
        {'eta':eta,'laurent_error':float(np.max(abs(
            second_order_outer_rate(f,.37,0.,eta,mu=.6,temperature=.32)-reconstruct_response(r,eta))))}
        for eta in (.08,.04,.02,.01)]
    data['paper_intrinsic_tau0']='NOT_VERIFIED'
    return data

if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True,indent=2))
