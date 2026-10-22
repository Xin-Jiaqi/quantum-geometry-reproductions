"""Round20 frozen-parameter hotspot, Fermi occupation and Kubo pole audit.

Local decomposition matches the existing cluster-interband Kubo protocol,
and compares cross-mu/same-side, near-resonance and small-gap pairs.
Net response is NEVER inferred from the absolute integrand alone.
"""
import json
from flux_lab.hotspot_diagnostics import weighted_k_audit,pair_resolved
from flux_lab.hexagonal_bz import hexagonal_symmetrized_model

def run():
    main=[
      weighted_k_audit("6",order=8,method="hexagon",components=((0,0),(0,1))),
      weighted_k_audit("6m'm'",order=8,method="hexagon",components=((0,0),(0,1))),
      weighted_k_audit("6",order=12,method="gauss",components=((0,0),(0,1)))
    ]
    H,V,S=hexagonal_symmetrized_model([.5145,-1.0456],"6")
    temps=[]
    for T in (.08,.17,.25):
        d=pair_resolved(H,V,S,mu=.13,temperature=T,z=.071+.004j)
        q=d["even"][0,0];same=d["components"]["same_side"][0,0]
        temps.append({"temperature":T,
                      "point_even_xx":[float(q.real),float(q.imag)],
                      "point_same_mu_side_xx":[float(same.real),float(same.imag)],
                      "min_intercluster_gap":d["min_intercluster_gap"]})
    return {"schema":"cui-flux-round20-kspace-hotspots-v1",
            "status":"model-only diagnostic; BZ convergence UNVERIFIED, full paper PARTIAL",
            "model":"synthetic triangular-periodic four-band orbital x physical spin",
            "main":main,"temperature_point_scan":temps,
            "important":"Forbidden tensor entries and absolute kernel are NOT converged observables; occupation labels overlap with resonance labels."}
if __name__=="__main__":print(json.dumps(run(),indent=2,ensure_ascii=False))
