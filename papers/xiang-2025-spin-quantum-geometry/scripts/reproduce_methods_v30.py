"""Deterministic finite-eta method JSON smoke, not intrinsic DC."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from methods.pt_doublet import hamiltonian,STAGGERED,projectors
from methods.finite_eta import second_order,resolvent
import numpy as np

if __name__=='__main__':
    k=(.38,.71);H,_=hamiltonian(*k)
    X=np.arange(16).reshape(4,4).astype(complex)
    r=second_order(lambda x,y:hamiltonian(x,y)[0],*k,STAGGERED,
                   mu=1.1,temperature=.17,eta=.1,step=.0004)
    res={'model':'synthetic_exact_PT_doublets','bands':np.linalg.eigvalsh(H).tolist(),
         'projector_ranks':[float(np.trace(P).real) for _,P in projectors(*k)],
         'finite_eta_value':r['value_real'],
         'solver_error':float(np.max(abs(resolvent(H,X,.1,1,0.,'band')-
                                       resolvent(H,X,.1,1,0.,'orbital')))),
         'protocol':r['protocol'],'scientific_scope':r['scientific_scope'],
         'full_paper_status':'PARTIAL'}
    print(json.dumps(res,sort_keys=True,indent=2))
