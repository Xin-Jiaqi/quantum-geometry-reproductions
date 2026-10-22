"""Audit exact Kramers-doublet equilibrium static Zeeman current Ward.

Independent oracles are a central difference in uniform B, a central
difference in k of equilibrium spin density, and spectral projectors.
This is an equilibrium static response; no finite-frequency KME claim.
"""
import json
import numpy as np
from flux_lab.degenerate_static_ward import (
    kramers_periodic_model,degenerate_static_ward,projector_static_ward,
    finite_zeeman_static,finite_k_spin_gradient)
from flux_lab.equilibrium_ward import static_zeeman_current_ward

PARAMS=dict(mu=.15,temperature=.22,kappa=.7)
def mx(x):return float(np.max(np.abs(x)))
def cmatrix(x):return {"real":np.asarray(x).real.tolist(),"imag":np.asarray(x).imag.tolist()}

def run():
    points=([0.,0.],[.31,-.22],[.57,.33])
    rows=[]
    for k in points:
        H,V,S=kramers_periodic_model(k)
        a=degenerate_static_ward(H,V,S,**PARAMS)
        b=projector_static_ward(H,V,S,**PARAMS)
        c=finite_zeeman_static(H,V,S,**PARAMS,step=1e-5)
        d=finite_k_spin_gradient(k,**PARAMS,step=1e-5)
        rows.append(dict(k=k,cluster_dimensions=a["cluster_dimensions"],
                         intra=cmatrix(a["intracluster"]),
                         inter=cmatrix(a["intercluster"]),
                         full=cmatrix(a["full_static"]),
                         projector_max_error=mx(a["full_static"]-b["full_static"]),
                         finite_B_max_error=mx(a["full_static"]-c),
                         spin_gradient_max_error=mx(a["full_static"]-d),
                         max_full=mx(a["full_static"]),
                         max_intra=mx(a["intracluster"])))
    return dict(schema="cui-round24-degenerate-static-ward-v1",
        model="triangular-periodic spin-orbital four-band with exact TR Kramers doublets at k=0",
        params=PARAMS,
        results=rows,
        scope="T>0 complete static equilibrium susceptibility; not finite-frequency metallic transport",
        original_paper="PARTIAL")

if __name__=="__main__":print(json.dumps(run(),indent=2))
