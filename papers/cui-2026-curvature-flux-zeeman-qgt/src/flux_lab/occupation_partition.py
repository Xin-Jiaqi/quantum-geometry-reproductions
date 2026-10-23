"""Pair-resolved k-space response diagnostics for Cui Flux synthetic models.

COMPLEMENT to the merged bz_hotspots module: band-pair occupations are
partitioned without altering the
response: cross-chemical-potential and same-side energy pairs, and subsets
near an external-frequency pole or a small intercluster gap. These subsets
overlap and must NOT all be summed together.

This intentionally distinguishes physical response magnitude from integration
cancellation. Near-degenerate gaps can have tiny f_n-f_m and negligible
observable weight. No metallic intraband, disorder, or SI interpretation.
"""
from __future__ import annotations
import numpy as np
from .nonabelian_response import _cluster_data
from .causal_response import _fermi
from .hexagonal_bz import hexagonal_bz_quadrature,hexagonal_symmetrized_model
from .primitive_bz import primitive_cell_quadrature
from .group_microscopic import build_group
from .core import PAULI

_PHYSICAL_SPIN=np.stack([np.kron(np.eye(2),s) for s in PAULI])

def pair_resolved(H,V,spin,*,mu=.13,temperature=.17,z=.071+.004j,
                  kappa=1.,gap_tolerance=1e-9,resonance_widths=2.,
                  small_gap_threshold=.02):
    """Pair-by-pair even/odd 2x3 tensors with explicit physical partitions.

    Cross-mu means (E_n-mu)*(E_m-mu)<0, NOT 'large Fermi difference'.
    Near-resonant means |+z-(E_m-E_n)|<=resonance_widths*|Im(z)| OR
    |-z-(E_m-E_n)|<=same threshold (both frequency branches).
    """
    if temperature<=0:raise ValueError("temperature>0 required for this audit")
    z=complex(z)
    if not np.isfinite(z.real) or not np.isfinite(z.imag) or z.imag<=0:
        raise ValueError("complex z must have positive imaginary regulator")
    if resonance_widths<=0 or small_gap_threshold<=0:
        raise ValueError("thresholds must be positive")
    d=_cluster_data(H,V,spin,gap_tolerance=gap_tolerance)
    e,v,s,labels=(d[k] for k in ("e","v","s","labels"))
    f=_fermi(e,mu,temperature)
    keys=("total","cross_mu","same_side","near_resonant","small_gap")
    parts={key:np.zeros((2,3),complex) for key in keys}
    pos=np.zeros((2,3),complex)
    neg=np.zeros((2,3),complex)
    min_gap=np.inf;min_detune=np.inf;resonant_count=0;resonant_max_df=0.
    transition_count=0
    for n in range(len(e)):
        for m in range(len(e)):
            if labels[n]==labels[m]:continue
            delta=e[m]-e[n]
            den_plus=z-delta
            den_minus=-z-delta
            if min(abs(den_plus),abs(den_minus))<1e-10:
                if abs(f[n]-f[m])>1e-10:raise ValueError("on-shell Kubo pole")
                continue
            numerator=-kappa*(f[n]-f[m])*np.outer(v[:,n,m],s[:,m,n])
            plus=numerator/den_plus;minus=numerator/den_minus
            value=(plus+minus)/2
            pos+=plus;neg+=minus
            parts["total"]+=value
            sector="cross_mu" if (e[n]-mu)*(e[m]-mu)<0 else "same_side"
            parts[sector]+=value
            detune=min(abs(den_plus),abs(den_minus))
            min_gap=min(min_gap,abs(delta))
            min_detune=min(min_detune,detune)
            transition_count+=1
            if detune<=resonance_widths*abs(z.imag):
                parts["near_resonant"]+=value
                resonant_count+=1
                resonant_max_df=max(resonant_max_df,abs(f[n]-f[m]))
            if abs(delta)<=small_gap_threshold:
                parts["small_gap"]+=value
    result=dict(even=parts["total"],odd=(pos-neg)/2,plus=pos,minus=neg,
                components=parts,energies=e,occupations=f,
                min_intercluster_gap=float(min_gap),
                min_pole_distance=float(min_detune),
                resonant_transition_count=resonant_count,
                max_occupation_difference_resonant=float(resonant_max_df),
                transition_count=transition_count)
    return result

def integrate_fermi_sectors(group,*,order=8,domain="hexagon",
                            method="gauss",component=(0,0),
                            mu=.13,temperature=.17,z=.071+.004j):
    """Exact thermal band-pair split integrated on an existing BZ quadrature.

    Not a second hotspot finder: use the accepted bz_hotspots module for
    node ranking, symmetry residuals and true spatial k-map visualization.
    This module answers a DIFFERENT question: which occupancy sectors
    survive *after* integration, including tiny Fermi-difference terms?
    """
    a,b=component
    if a not in (0,1) or b not in (0,1,2):raise ValueError("component=(0/1,0/1/2)")
    ops=build_group(group)
    if domain=="hexagon":
        k,w,area=hexagonal_bz_quadrature(order)
    elif domain=="primitive":
        kw={"shift":(.23,-.17)} if method=="midpoint" else {}
        k,w,area=primitive_cell_quadrature(order,method=method,**kw)
    else:raise ValueError("domain must be hexagon or primitive")
    names=("total","cross_mu","same_side","near_resonant","small_gap")
    totals={name:np.zeros((2,3),complex) for name in names}
    absavg={name:np.zeros((2,3),float) for name in names}
    min_gap=np.inf;res_max_df=0.;n_res=0
    for point,wt in zip(k,w):
        H,V,S=hexagonal_symmetrized_model(point,ops)
        d=pair_resolved(H,V,S,mu=mu,temperature=temperature,z=z)
        for name in names:
            x=d["components"][name]
            totals[name]+=wt*x
            absavg[name]+=wt*abs(x)
        min_gap=min(min_gap,d["min_intercluster_gap"])
        res_max_df=max(res_max_df,d["max_occupation_difference_resonant"])
        n_res+=d["resonant_transition_count"]
    v=totals["total"][a,b];mod=absavg["total"][a,b]
    def complex_pair(x):
        return [float(x.real),float(x.imag)]
    sectors={}
    for name in names:
        sectors[name]={"net":complex_pair(totals[name][a,b]),
                       "net_abs":float(abs(totals[name][a,b])),
                       "weighted_absolute":float(absavg[name][a,b])}
    return {"group":group,"order":int(order),"domain":domain,"method":method,
            "component":[int(a),int(b)],"nodes":len(w),"bz_area":float(area),
            "params":{"mu":mu,"temperature":temperature,
                      "z":[float(z.real),float(z.imag)]},
            "sectors":sectors,
            "cancellation_ratio":float(mod/max(abs(v),1e-30)),
            "minimum_sampled_intercluster_gap":float(min_gap),
            "max_df_of_near_pole_pairs":float(res_max_df),
            "near_pole_pair_visits":int(n_res),
            "complete_partition_residual":float(np.max(abs(
                totals["total"]-totals["cross_mu"]-totals["same_side"]))),
            "scope":"thermal band-pair sector diagnostic only; no BZ convergence claim"}
