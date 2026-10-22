"""Round23: equilibrium no-current Ward test, independent finite-B oracle."""
import json
import numpy as np
from flux_lab.core import PAULI
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model
from flux_lab.equilibrium_ward import (static_zeeman_current_ward,
                                      zeeman_finite_difference,
                                      spin_k_gradient_fd,static_grid_audit)

def cplx(a):
    a=np.asarray(a)
    return {"real":a.real.tolist(),"imag":a.imag.tolist()}

def main():
    k=np.array([.47,-.31]);model=lambda k:hexagonal_symmetrized_model(k,"1")
    H,V,S=model(k)
    r=static_zeeman_current_ward(H,V,S,mu=2.3,temperature=.3)
    b=zeeman_finite_difference(H,V,S,mu=2.3,temperature=.3)
    dk=spin_k_gradient_fd(model,k,mu=2.3,temperature=.3)
    grids={}
    for n in (4,8,12):
        x=static_grid_audit("1",order=n,mu=2.3,temperature=.3)
        grids[str(n)]={"interband":cplx(x["interband"]),
                       "occupation_derivative":cplx(x["occupation_derivative"]),
                       "static_total":cplx(x["full_static"]),
                       "max_static_total":float(np.max(abs(x["full_static"]))),
                       "exact_discrete_periodic_divergence_max":
                          float(np.max(abs(x["telescoping_derivative"])))}
    return {"schema":"cui-flux-round23-equilibrium-ward-v1",
            "status":"equilibrium fixed-basis static Ward, not finite-frequency KME",
            "point_k":k.tolist(),"mu":2.3,"temperature":.3,
            "point_interband":cplx(r["interband"]),
            "point_occupation_derivative":cplx(r["occupation_derivative"]),
            "point_total":cplx(r["full_static"]),
            "point_finiteB_error":float(np.max(abs(r["full_static"]-b))),
            "point_spin_gradient_error":float(np.max(abs(r["full_static"]-dk))),
            "primitive_midpoint_grid":grids,
            "physical_claim":"exact equilibrium total response is a periodic Bloch spin-density derivative; its continuum BZ integral is zero. Discrete spectral quadrature requires convergence.",
            "full_paper":"PARTIAL"}

if __name__=="__main__":print(json.dumps(main(),indent=2))
