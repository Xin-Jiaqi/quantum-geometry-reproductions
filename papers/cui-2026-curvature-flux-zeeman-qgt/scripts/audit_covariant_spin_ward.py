"""Reproducible Round25 k-dependent spin/Bloch-gauge static Ward audit."""
import json
import numpy as np
from flux_lab.covariant_spin_ward import (
 periodic_spin_embedding,covariant_static_ward,finite_zeeman_current_derivative,
 finite_k_spin_density_gradient)

def maximum(A):return float(np.max(abs(A)))
def run():
    rows=[]
    for k in ([0.,0.],[.36,-.29],[.71,.23]):
        row={"k":list(k)}
        for v in (0.,.16):
            H,V,S,dS,A=periodic_spin_embedding(k,effective_spin_variation=v,rotate=True)
            out=covariant_static_ward(H,V,S,dS,A)
            finiteB=finite_zeeman_current_derivative(H,V,S,dS,A)
            finiteK=finite_k_spin_density_gradient(np.asarray(k),effective_spin_variation=v)
            H0,V0,S0,dS0,A0=periodic_spin_embedding(k,effective_spin_variation=v,rotate=False)
            ref=covariant_static_ward(H0,V0,S0,dS0,A0)
            row[f"variation_{v:g}"]={
               "full_static_real":out["full_static"].real.tolist(),
               "full_static_imag":out["full_static"].imag.tolist(),
               "gauge_invariance_max_error":maximum(out["full_static"]-ref["full_static"]),
               "finite_B_error":maximum(out["full_static"]-finiteB),
               "finite_k_error":maximum(out["full_static"]-finiteK),
               "explicit_spin_contact_max":maximum(out["spin_vertex_contact"]),
               "population_only_error":maximum(out["occupation"]-finiteK),
               "raw_vs_covariant_velocity_max":maximum(V-out["covariant_velocity"])
            }
        rows.append(row)
    return {"schema":"cui-flux-round25-kdependent-spin-covariant-ward-v1",
            "physics":"static finite-T equilibrium only; full paper PARTIAL",
            "conditions":{"mu":.15,"temperature":.22,"kappa":.7},
            "rows":rows,
            "scope":"synthetic Bloch gauge, synthetic effective spin variation, no material/SI/dynamical KME"}
if __name__=="__main__":print(json.dumps(run(),indent=2,ensure_ascii=False))
