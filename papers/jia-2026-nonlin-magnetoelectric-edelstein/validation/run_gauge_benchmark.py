"""Benchmark intentionally uses a SYNTHETIC 4-orbital insulator, not a material."""
import sys,csv,json
from pathlib import Path
from tests.test_covariant import toy, physical_model, twist
from nmee_lab import mesh,evaluate_insulating_static
import numpy as np
root=Path(__file__).resolve().parent
TB=toy();m,spin,A=physical_model(TB,True);mz,_,_=physical_model(TB,False)
mt=twist(m,np.diag([0.,.4,0.,.4]))
rows=[]
for n in (11,21,41,81):
    k=mesh(TB.reciprocal_nm,n,offset=(.31,.27))
    for B in (.1,.05,.025):
        kw=dict(observable='top_z',magnetic='x',beta=0,n_occupied=2,B_step_T=B)
        x=evaluate_insulating_static(m,k,**kw)
        y=evaluate_insulating_static(mt,k,**kw)
        z=evaluate_insulating_static(mz,k,**kw)
        rows.append(dict(n=n,B_T=B,gamma_nonzero_A=x['response'],gamma_k_gauge_twist=y['response'],gamma_A0=z['response'],
                         gauge_abs_err=abs(x['response']-y['response']),connection_effect=x['response']-z['response'],min_gap_eV=x['min_gap_eV']))
with (root/'synthetic_covariant_benchmark.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
meta=dict(min_gap_eV=min(r['min_gap_eV'] for r in rows),max_gauge_abs_err=max(r['gauge_abs_err'] for r in rows),
          synthetic_only=True,meshes=[11,21,41,81],field_steps_T=[.1,.05,.025],
          response_definition='T=0 insulator; occupied-unoccupied static interband chi derivative d/dB, with D_k=partial_k-i[A_k, .]')
(root/'synthetic_benchmark_meta.json').write_text(json.dumps(meta,indent=2))
print(meta)
print('last n=81 B=.025',rows[-1])
