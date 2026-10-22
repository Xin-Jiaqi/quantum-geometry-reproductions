"""Pair-resolved k-space response diagnostics for Cui Flux synthetic models.

Pointwise clean cluster-interband Kubo is partitioned *without altering* the
response: cross-chemical-potential and same-side energy pairs, and subsets
near an external-frequency pole or a small intercluster gap. These subsets
overlap and must NOT all be summed together.

This intentionally distinguishes physical response magnitude from integration
cancellation. Near-degenerate gaps can have tiny f_n-f_m and negligible
observable weight. No metallic intraband, disorder, or SI interpretation.
"""
from __future__ import annotations
import numpy as np
from .nonabelian_response import _cluster_data,cluster_interband_kubo
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
    odd=np.zeros((2,3),complex)
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

def weighted_k_audit(group,*,order=8,method="hexagon",components=((0,0),(0,1)),
                     mu=.13,temperature=.17,z=.071+.004j,
                     resonance_widths=2.,small_gap_threshold=.02,top_n=5):
    """Integration and hotspot diagnostics, with no symmetry projection.

    Intrinsic absolute-kernel integral = sum_i w_i |chi_even_ab(k_i)|,
    cancellation ratio = that integral / |sum_i w_i chi_even_ab(k_i)|.
    As net tends to zero the ratio diverges and is NOT a condition number
    certified against machine precision.
    """
    ops=build_group(group)
    if method=="hexagon":
        pts,w,area=hexagonal_bz_quadrature(order)
    elif method in ("gauss","midpoint"):
        kw={"shift":(.23,-.17)} if method=="midpoint" else {}
        pts,w,area=primitive_cell_quadrature(order,method=method,**kw)
    else:
        raise ValueError("method must be hexagon, gauss or midpoint")
    comps=[]
    for a,b in components:
        if not (0<=a<2 and 0<=b<3):
            raise ValueError("spin-current index requires 2x3 response")
        comps.append((int(a),int(b)))
    names=("total","cross_mu","same_side","near_resonant","small_gap")
    totals={name:np.zeros((2,3),complex) for name in names}
    absolutes={name:np.zeros((2,3),float) for name in names}
    min_gap=np.inf;min_pole=np.inf;max_resonant_df=0.;num_resonant=0
    sampled=[]
    for point,weight in zip(pts,w):
        H,V=hexagonal_symmetrized_model(point,ops)[:2]
        d=pair_resolved(H,V,_PHYSICAL_SPIN,
                 mu=mu,temperature=temperature,z=z,
                 resonance_widths=resonance_widths,
                 small_gap_threshold=small_gap_threshold)
        for name in names:
            x=d["components"][name]
            totals[name]+=weight*x
            absolutes[name]+=weight*abs(x)
        min_gap=min(min_gap,d["min_intercluster_gap"])
        min_pole=min(min_pole,d["min_pole_distance"])
        num_resonant+=d["resonant_transition_count"]
        max_resonant_df=max(max_resonant_df,d["max_occupation_difference_resonant"])
        sampled.append((np.asarray(point).copy(),float(weight),d))
    output={}
    for a,b in comps:
        net=totals["total"][a,b];mag=absolutes["total"][a,b]
        top=sorted(sampled,key=lambda row:row[1]*abs(row[2]["even"][a,b]),
                   reverse=True)[:top_n]
        output[f"{a}{b}"]={
            "net":[float(net.real),float(net.imag)],
            "net_abs":float(abs(net)),
            "absolute_integral":float(mag),
            "cancellation_ratio":float(mag/max(abs(net),1e-30)),
            "cross_mu_net":[float(totals["cross_mu"][a,b].real),
                            float(totals["cross_mu"][a,b].imag)],
            "same_side_net":[float(totals["same_side"][a,b].real),
                             float(totals["same_side"][a,b].imag)],
            "cross_mu_absolute_integral":float(absolutes["cross_mu"][a,b]),
            "same_side_absolute_integral":float(absolutes["same_side"][a,b]),
            "near_resonant_absolute_integral":float(absolutes["near_resonant"][a,b]),
            "small_gap_absolute_integral":float(absolutes["small_gap"][a,b]),
            "top_nodes":[{"k":p.tolist(),"quadrature_weight":wt,
                          "weighted_abs_kernel":float(wt*abs(d["even"][a,b])),
                          "min_gap":d["min_intercluster_gap"],
                          "min_pole_distance":d["min_pole_distance"]}
                         for p,wt,d in top]}
    return {
        "group":group,"order":int(order),"method":method,"nodes":len(pts),
        "bz_area":float(area),"weight_sum":float(np.sum(w)),
        "parameters":{"mu":mu,"temperature":temperature,
                      "frequency":[float(z.real),float(z.imag)],
                      "resonance_eta_multiples":resonance_widths,
                      "small_gap_threshold":small_gap_threshold},
        "min_intercluster_gap":float(min_gap),
        "min_pole_distance":float(min_pole),
        "max_resonant_occupation_difference":float(max_resonant_df),
        "resonant_transition_visits":int(num_resonant),
        "max_total_partition_residual":float(np.max(abs(totals["total"]-
                 totals["cross_mu"]-totals["same_side"]))),
        "components":output,
        "status":"diagnostic only; no finite-order BZ convergence claim"}
