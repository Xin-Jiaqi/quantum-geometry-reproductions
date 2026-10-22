"""Round18 convergence audit for hexagonal first-BZ response quadrature.

Orders are a convergence study on the synthetic periodic model only.
The highest order is a numerical reference, not an exact analytic answer.
"""
import json
import numpy as np
from flux_lab.hexagonal_bz import hexagonal_bz_response
from flux_lab.group_microscopic import project_frequency_tensor

ORDERS=(2,3,4,5,6,8,10,12,16,20,24)
GROUPS=("6","6m'm'")

def maxabs(x):
    return float(np.max(abs(x)))

def serial_matrix(x):
    return {"real":np.asarray(x).real.tolist(),"imag":np.asarray(x).imag.tolist()}

def run():
    output={}
    for group in GROUPS:
        results={}
        for order in ORDERS:
            d=hexagonal_bz_response(group,order=order,z=.071+.004j,mu=.13,temperature=.17)
            results[order]=d
        ref=results[ORDERS[-1]]
        rows=[]
        for order in ORDERS:
            d=results[order]
            rows.append({
                "order":order,
                "nodes":len(d["weights"]),
                "max_even":maxabs(d["even"]),
                "max_odd":maxabs(d["odd"]),
                "even_delta_vs_order24":maxabs(d["even"]-ref["even"]),
                "odd_delta_vs_order24":maxabs(d["odd"]-ref["odd"]),
                "even_projector_residual":maxabs(d["even"]-project_frequency_tensor(group,d["even"],"even")),
                "odd_projector_residual":maxabs(d["odd"]-project_frequency_tensor(group,d["odd"],"odd")),
                "even":serial_matrix(d["even"]),
                "odd":serial_matrix(d["odd"])
            })
        output[group]=rows
    return {
        "schema":"cui-flux-round18-hexagonal-bz-convergence-v3",
        "model":"synthetic periodic triangular reciprocal model",
        "frequency":[.071,.004],
        "orders":list(ORDERS),
        "interpretation":"order-24 is a numerical reference only; no real-material or SI claim",
        "groups":output
    }

if __name__=="__main__":
    print(json.dumps(run(),indent=2,ensure_ascii=False))
