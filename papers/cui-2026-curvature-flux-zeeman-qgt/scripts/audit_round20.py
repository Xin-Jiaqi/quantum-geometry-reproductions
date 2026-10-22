"""Round20: k-resolved spectral gap, occupancy, and cancellation audit.

A large unprojected primitive-cell tensor error can be either high local
spectral variation or cancellation. Measure both; do not declare a mechanism
from the group projector residual or a few finite quadratures. All
frequencies/temperatures are model natural units and synthetic only.
"""
from __future__ import annotations
import argparse,json
import numpy as np
from flux_lab.bz_hotspots import bz_hotspot_audit

GROUPS=("6","6m'm'")
# Fixed Cartesian observables: C6 diagonal yy; C6 antiunitary xy.
# Avoid sampling-dependent or machine-dependent auto argmax among symmetry
# equivalent tensor components when interpreting cancellation/hotspots.
COMPONENTS={"6":(1,1),"6m'm'":(0,1)}

def compact(d):
    A=d["even"];B=d["odd"]
    def pack(C):
        return dict(net={"real":C["net"].real.tolist(),"imag":C["net"].imag.tolist()},
            max_abs_net=float(np.max(abs(C["net"]))),
            abs_integral_selected=C["selected_abs_integral"],
            retention=C["selected_retention"],
            cancellation=C["selected_cancellation"],
            concentration_top5pct=C["top_share"],
            max_local_abs=C["maximum_local_abs"])
    return dict(domain=d["domain"],method=d["method"],order=d["order"],
        nodes=d["nodes"],component=list(d["component"]),mu=d["mu"],
        temperature=d["temperature"],z=[d["z"].real,d["z"].imag],
        min_intercluster_gap=d["min_gap"],
        min_active_resonance_distance=d["min_active_resonance_distance"],
        max_pair_occupation_difference=d["max_occupation_difference"],
        active_directed_pair_evaluations=d["total_active_directed_transition_count"],
        even=pack(A),odd=pack(B),hotspots=d["hotspots"])


def run(*,deep=False):
    group_rows={}
    for group in GROUPS:
        hex_order=6 if deep else 5
        component=COMPONENTS[group]
        anchor=bz_hotspot_audit(group,domain="hexagon",order=hex_order,
                                component=component,hotspot_count=6)
        other_hex=bz_hotspot_audit(group,domain="hexagon",order=4 if deep else 3,
                                    component=component,hotspot_count=6)
        prim_gauss=bz_hotspot_audit(group,domain="primitive",method="gauss",
                                  order=16 if deep else 10,
                                  component=component,hotspot_count=6)
        prim_mid=bz_hotspot_audit(group,domain="primitive",method="midpoint",
                                 shift=(.23,-.17),order=20 if deep else 12,
                                 component=component,hotspot_count=6)
        samples=[compact(v) for v in [other_hex,anchor,prim_gauss,prim_mid]]
        for s in samples:
            even=np.array(s["even"]["net"]["real"])+1j*np.array(s["even"]["net"]["imag"])
            s["even_full_tensor_delta_vs_hex_anchor"]=float(np.max(abs(even-anchor["even"]["net"])))
            odd=np.array(s["odd"]["net"]["real"])+1j*np.array(s["odd"]["net"]["imag"])
            s["odd_full_tensor_delta_vs_hex_anchor"]=float(np.max(abs(odd-anchor["odd"]["net"])))
        group_rows[group]={"selected_component":list(component),
                           "anchor":"hexagon order "+str(hex_order),
                           "samples":samples}
    temp=[]
    for T in ((.1,.17,.35,.6) if deep else (.1,.17,.35)):
        item=bz_hotspot_audit("6",domain="hexagon",order=4,temperature=T,
                              component=tuple(group_rows["6"]["selected_component"]))
        temp.append(compact(item))
    return dict(schema="cui-flux-round20-local-gap-hotspot-cancellation-v1",
        groups=group_rows,temperature_scan_group6=temp,
        scientific_scope="synthetic clean cluster interband only; no material or SI",
        warnings=[
          "Large cancellation fraction indicates net response small relative to integral of modulus; not a proof of singular Berry geometry",
          "Minimum intercluster gap and resonance distance are discrete sampled values, not global extrema",
          "Hotspots depend on sampling and quadrature weights; convergence still requires independent refined grids",
          "All tensor components retained unprojected; no automatic convergence PASS"
        ])

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--deep",action="store_true")
    a=p.parse_args()
    print(json.dumps(run(deep=a.deep),ensure_ascii=False,indent=2))
