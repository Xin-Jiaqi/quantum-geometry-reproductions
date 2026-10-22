"""Round22 cancellation conditioning and empirical multi-grid BZ uncertainty.

WARNING: independence across finite grids is useful but NOT a theorem that
their observed spread upper-bounds the unknown integration error.
The only mathematical inequality here is the explicitly conditional
amplification formula using an assumed weighted L1 integrand error.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from flux_lab.bz_error_budget import (cancellation_condition,
                                      shifted_midpoint_ensemble)

ROOT=Path(__file__).resolve().parents[1]
MODEL_CHANNELS={"6":(1,1),"6m'm'":(0,1)}

def complex_data(z):
    z=complex(z)
    return [float(z.real),float(z.imag)]

def summarize(group,diag):
    sample=next(x for x in diag["groups"][group]["samples"]
                if x["domain"]=="hexagon" and x["order"]==5)
    retention=sample["even"]["retention"]
    return {"sample_order":5,"observed_retention":retention,
            "required_for_1_percent_relative_net":
            cancellation_condition(retention,desired_relative_net_error=.01)}

def run(deep=False):
    frozen=json.loads((ROOT/"results/frozen/round20_khotspot_cancellation.json").read_text())
    settings={"coarse_order":12,"fine_order":18,"hex_order":8,"gauss_order":18} if deep else {
        "coarse_order":6,"fine_order":9,"hex_order":5,"gauss_order":9}
    groups={}
    for group,comp in MODEL_CHANNELS.items():
        d=shifted_midpoint_ensemble(group,component=comp,parity="even",**settings)
        c=d["comparison"]
        groups[group]={
          "observed_component":list(comp),
          "symmetry_allowed":d["symmetry_allowed"],
          "sampling_orders":settings,
          "num_shifted_midpoint_grids":len(d["fine_midpoint"]),
          "shifted_midpoint_coarse":[complex_data(x) for x in d["coarse_midpoint"]],
          "shifted_midpoint_fine":[complex_data(x) for x in d["fine_midpoint"]],
          "finite_hex_value":complex_data(c["finite_hex_anchor"]),
          "finite_primitive_gauss":complex_data(c["finest_primitive_gauss"]),
          "fine_midpoint_mean":complex_data(c["finest_midpoint_shift_mean"]),
          "within_fine_shift_absolute_spread":c["fine_shift_spread"],
          "mean_midpoint_refinement_change":c["midpoint_refinement_difference"],
          "across_finest_methods_spread":c["independent_method_spread"],
          "empirical_absolute_disagreement":c["empirical_disagreement"],
          "empirical_relative_disagreement_to_finite_hex":
                c["empirical_disagreement_relative_to_finite_anchor"],
          "hex_projector_residual":d["hex_projector_residual"],
          "primitive_gauss_projector_residual":d["gauss_projector_residual"],
          "convergence_status":c["numerical_convergence_status"],
          "mathematical_condition_from_prior_round20":summarize(group,frozen),
        }
    return dict(
      schema="cui-flux-round22-cancellation-conditioned-bz-uncertainty-v1",
      type="synthetic clean cluster-interband model only",
      numerical_convergence="UNVERIFIED",
      mathematical_result="Conditional: if integral |delta F| <= eps integral |F|, relative net error <= eps / retention",
      empirical_result="Finite-sampling spread is NOT a rigorous or statistical confidence interval",
      method="Cartesian fixed component, four independent shifted midpoint meshes at two resolutions, primitive Gauss and D6 hexagon Gauss",
      groups=groups,
      negative_control="Periodic Fourier aliasing at n=grid size defeats an unshifted midpoint rule and is checked by unit tests",
      full_paper="PARTIAL")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--deep",action="store_true")
    a=p.parse_args();print(json.dumps(run(a.deep),indent=2,ensure_ascii=False))
