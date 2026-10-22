"""Round19: independent triangular primitive-cell versus hexagonal Wigner-Seitz audit.

Compute *unprojected* reciprocal-parallelogram Gauss and shifted MP grids
and compare their full complex even/odd tensors against the separately
triangulated, D6-symmetrized first-hexagonal-BZ quadrature.

The highest order from any sequence is a NUMERICAL COMPARATOR, not an exact
physical answer. Symmetry projector residual != quadrature convergence.
No real Wannier, collision or SI claim.
"""
from __future__ import annotations
import argparse,json
import numpy as np
from flux_lab.hexagonal_bz import hexagonal_bz_response
from flux_lab.primitive_bz import primitive_bz_response
from flux_lab.group_microscopic import project_frequency_tensor

GROUPS=("6","6m'm'")
def maxabs(a):
    return float(np.max(abs(a)))

def encode(arr):
    return {"real":np.asarray(arr).real.tolist(),"imag":np.asarray(arr).imag.tolist()}

def one_result(group,kind,order,ref):
    if kind=="hexagon":
        d=hexagonal_bz_response(group,order=order)
        N=len(d["weights"])
    else:
        opts={"method":kind}
        if kind=="midpoint":opts["shift"]=(.23,-.17)
        d=primitive_bz_response(group,order=order,**opts)
        N=d["number_points"]
    ev,od=d["even"],d["odd"]
    return dict(method=kind,order=order,nodes=N,
        even_max=maxabs(ev),odd_max=maxabs(od),
        even_diff_vs_hex_ref=maxabs(ev-ref["even"]),
        odd_diff_vs_hex_ref=maxabs(od-ref["odd"]),
        even_projector_residual=maxabs(ev-project_frequency_tensor(group,ev,"even")),
        odd_projector_residual=maxabs(od-project_frequency_tensor(group,od,"odd")),
        even=encode(ev),odd=encode(od))

def run(deep=False):
    hex_orders=(8,12,16) if deep else (8,12)
    gauss_orders=(8,12,20,32) if deep else (8,12,20)
    mp_orders=(12,20,32,48) if deep else (12,20,32)
    data={}
    for group in GROUPS:
        hex_ref=hexagonal_bz_response(group,order=hex_orders[-1])
        rows=[]
        for method,orders in [("hexagon",hex_orders),("gauss",gauss_orders),
                              ("midpoint",mp_orders)]:
            for order in orders:
                rows.append(one_result(group,method,order,hex_ref))
        data[group]={"reference_hex_order":hex_orders[-1],
                     "reference_hex_nodes":len(hex_ref["weights"]),
                     "rows":rows}
    return {"schema":"cui-flux-round19-independent-primitive-bz-v1",
            "scope":"independent periodic primitive-cell integration vs first hexagon",
            "parameters":{"z":[.071,.004],"mu":.13,"temperature":.17,
                          "reciprocal_cell":"B times [-pi,pi)^2",
                          "midpoint_shift":[.23,-.17]},
            "deep":bool(deep),"groups":data,
            "interpretation":"No automatic convergence PASS; compare methods and full tensor components. The highest-order hex value is finite-order only."}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--deep",action="store_true")
    a=ap.parse_args()
    print(json.dumps(run(a.deep),indent=2,ensure_ascii=False))
