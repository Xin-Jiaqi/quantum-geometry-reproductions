"""Reproducible Round17 non-Abelian + hexagonal-BZ acceptance summary."""
import json
import numpy as np
from flux_lab.core import PAULI
from flux_lab.nonabelian_response import (
    cluster_interband_kubo,cluster_interband_liouville,cluster_zeeman_qgt,cluster_structure)
from flux_lab.hexagonal_bz import (
    hexagonal_bz_quadrature,hexagonal_local_residual,hexagonal_bz_response)
from flux_lab.group_microscopic import project_frequency_tensor

I=np.eye(2,dtype=complex)
SPIN=np.stack([np.kron(I,s) for s in PAULI])

def kramers_fixture():
    H=1.2*np.kron(PAULI[2],I)+.31*np.kron(PAULI[1],PAULI[0])
    V=np.stack([
        .55*np.kron(PAULI[0],PAULI[0])+.21*np.kron(PAULI[1],I),
        .47*np.kron(PAULI[0],PAULI[1])+.13*np.kron(PAULI[0],PAULI[2])])
    return H,V,SPIN

def maximum(*xs):
    return float(max(np.max(abs(x)) for x in xs))

def run():
    H,V,S=kramers_fixture();z=.071+.009j
    a=cluster_interband_kubo(H,V,S,z=z,mu=0.,temperature=.17)
    b=cluster_interband_liouville(H,V,S,z=z,mu=0.,temperature=.17)
    q=cluster_zeeman_qgt(H,V,S,cluster_index=0)
    points,weights,area=hexagonal_bz_quadrature(4)
    local={g:hexagonal_local_residual(g) for g in ("6","6mm","6m'm'","3m'")}
    bz={}
    for g in ("6","6m'm'"):
        d=hexagonal_bz_response(g,order=3)
        bz[g]={
            "even_projector_residual":float(np.max(abs(d["even"]-project_frequency_tensor(g,d["even"],"even")))),
            "odd_projector_residual":float(np.max(abs(d["odd"]-project_frequency_tensor(g,d["odd"],"odd")))),
            "max_even":float(np.max(abs(d["even"]))),
            "max_odd":float(np.max(abs(d["odd"]))),
            "order":d["order"],"points":len(d["weights"])}
    return {
        "schema":"cui-flux-round17-nonabelian-hex-bz-v1",
        "scope":"synthetic model-method acceptance only; full paper PARTIAL",
        "nonabelian":{
            "cluster_dimensions":cluster_structure(H)["dimensions"],
            "kubo_liouville_max_residual":maximum(a[0]-b[0],a[1]-b[1]),
            "max_qgt_trace":float(np.max(abs(q["trace"]))),
            "drive_frequency":[z.real,z.imag],
            "protocol":"clean interband; near-degenerate clusters flattened to their mean energy"},
        "hexagonal_bz":{
            "reciprocal_basis": [[1.,.5],[0.,float(np.sqrt(3.)/2)]],
            "first_bz_area":area,"quadrature_order":4,
            "quadrature_points":len(points),"weight_sum":float(weights.sum()),
            "local_covariance_maxima":{g:maximum(*[np.array(v) for k,v in r.items() if k in ("H","velocity","chi","alpha")]) for g,r in local.items()},
            "integrated_projector_checks":bz,
            "quadrature":"six triangles averaged with mirror image to enforce full D6 node symmetry"},
        "limitations":["not a real material","not a conserving intraband/metallic response",
            "no finite scattering, contacts, SI conversion, or Table I adjudication",
            "order-3 BZ integrations are symmetry checks, not convergence claims"]
    }

if __name__=="__main__":
    print(json.dumps(run(),indent=2,ensure_ascii=False))
