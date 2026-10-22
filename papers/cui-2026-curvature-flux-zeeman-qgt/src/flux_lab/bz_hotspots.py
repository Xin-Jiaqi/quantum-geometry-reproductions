"""Diagnose cancellation, spectral gaps and local k hotspots in Cui Flux BZ.

The object being measured is the synthetic model's clean cluster-*interband*
JB tensor, never a physical intraband metal conductivity. All node sampling
schemes are existing independently validated quadratures; no group projector
is applied. Symmetry-forbidden components can have zero net by construction
but large absolute local response. Report BOTH complex net and L1 magnitude.

For delta = e_m-e_n, transition q = f_n-f_m, A=q*v_nm*s_mn,
chi_even = -kappa sum A*delta/(z**2-delta**2)
chi_odd  = -kappa sum A*z/(z**2-delta**2).
This is a DIFFERENT algebraic evaluation than separately computing chi(z),
chi(-z) and halving. Any comparison to two calls of the original Kubo is an
independent numerical oracle, not an additional fitted identity.
"""
from __future__ import annotations
import numpy as np
from .nonabelian_response import _cluster_data
from .causal_response import _fermi
from .hexagonal_bz import (INV_BASIS,hexagonal_bz_quadrature,
                           hexagonal_symmetrized_model)
from .primitive_bz import primitive_cell_quadrature
from .group_microscopic import build_group
from .hexagonal_bz import hexagonal_compatible


def local_transition_diagnostic(H,V,spin,*,mu=.13,temperature=.17,z=.071+.004j,
                                kappa=1.,gap_tolerance=1e-9,occupation_cutoff=1e-5):
    """One diagonalization; directed band pairs and diagnostic denominators.

    The 'resonance_distance' diagnostic is the min |abs(delta)-abs(Re z)|
    among cluster pairs weighted by a nonzero occupation difference; it is
    NOT a linewidth-corrected spectral function, and does not substitute
    for a true real-frequency Kubo calculation.
    """
    D=_cluster_data(H,V,spin,gap_tolerance=gap_tolerance)
    e=D["e"];labels=D["labels"];f=_fermi(e,mu,temperature);z=complex(z)
    if temperature<0 or not np.isfinite([z.real,z.imag]).all():
        raise ValueError("finite drive and T>=0 required")
    if not np.isfinite(occupation_cutoff) or occupation_cutoff<0:
        raise ValueError("occupation_cutoff must be finite and nonnegative")
    even=np.zeros((2,3),complex);odd=np.zeros((2,3),complex)
    gaps=[];active=[];occ=[];transitions=0
    for n in range(len(e)):
        for m in range(len(e)):
            if labels[n]==labels[m]:continue
            de=float(e[m]-e[n]);gaps.append(abs(de))
            df=float(f[n]-f[m]);occ.append(abs(df))
            if abs(df)<=occupation_cutoff:continue
            active.append(abs(abs(de)-abs(z.real)))
            denom=z*z-de*de
            if abs(denom)<1e-12:raise ValueError("on-shell pole; use nonzero eta")
            amp=(-kappa*df)*np.outer(D["v"][:,n,m],D["s"][:,m,n])/denom
            even+=amp*de;odd+=amp*z
            transitions+=1
    if not gaps:raise ValueError("Only one cluster: no interband transitions")
    return dict(even=even,odd=odd,
                min_intercluster_gap=float(min(gaps)),
                closest_active_resonance_distance=float(min(active)) if active else None,
                active_transition_count=transitions,
                max_pair_occupation_difference=float(max(occ)),
                band_energies=[float(x) for x in e],
                cluster_dimensions=[int(np.count_nonzero(labels==n))
                                    for n in range(int(labels[-1])+1)])


def field_diagnostics(values,weights,*,component=None,top_fraction=.05):
    """Compute net, integrated modulus, and k-space cancellation for a tensor.

    values: (Nk,2,3) complex, normalized positive weights: (Nk,)
    component=None automatically chooses a,b with largest |integrated|.
    top_share reports share of weighted |local integrand| in the top
    ceil(top_fraction*Nk) sampled nodes; unequal quadrature weights matter.
    """
    val=np.asarray(values,complex);w=np.asarray(weights,float)
    if val.ndim!=3 or val.shape[1:]!=(2,3) or w.shape!=(len(val),):
        raise ValueError("expected (Nk,2,3) tensor field and (Nk,) weights")
    if len(w)==0 or not np.all(np.isfinite(val)) or not np.all(np.isfinite(w)):
        raise ValueError("nonfinite/empty field")
    if np.any(w<=0) or not np.isclose(w.sum(),1.,atol=1e-10):
        raise ValueError("weights must be positive and normalized")
    if not 0<top_fraction<=1:raise ValueError("top_fraction must be in (0,1]")
    net=np.einsum("n,nab->ab",w,val,optimize=True)
    l1=np.einsum("n,nab->ab",w,np.abs(val),optimize=True)
    if component is None:
        ab=np.unravel_index(np.argmax(abs(net)),net.shape)
    else:
        if len(component)!=2 or component[0] not in (0,1) or component[1] not in (0,1,2):
            raise ValueError("component must be (a=0/1,b=0/1/2)")
        ab=tuple(component)
    contribution=w*np.abs(val[:,ab[0],ab[1]])
    order=np.argsort(contribution)[::-1]
    n_top=max(1,int(np.ceil(top_fraction*len(w))))
    L=float(l1[ab]);N=float(abs(net[ab]))
    return dict(net=net,integrated_abs=l1,component=ab,
                selected_net=complex(net[ab]),selected_abs_integral=L,
                selected_retention=N/L if L>0 else 1.,
                selected_cancellation=1.-N/L if L>0 else 0.,
                top_fraction=float(top_fraction),top_share=float(np.sum(contribution[order[:n_top]])/L)
                       if L>0 else 0.,
                top_indices=[int(i) for i in order[:n_top]],
                maximum_local_abs=float(np.max(np.abs(val[:,ab[0],ab[1]]))))


def bz_hotspot_audit(group,*,domain="hexagon",order=6,method="gauss",
                     shift=(0.,0.),mu=.13,temperature=.17,z=.071+.004j,
                     component=None,hotspot_count=6,top_fraction=.05):
    """Quantify BZ cancellation and spectral gaps without enforcing symmetry.

    Returns one response tensor per frequency parity and the dominant
    hotspots for the selected even-frequency component. Sampling weights
    and k-space positions are included to diagnose nonuniform quadratures.
    """
    if not isinstance(hotspot_count,(int,np.integer)) or not 1<=hotspot_count<=30:
        raise ValueError("hotspot_count must be integer in 1..30")
    ops=build_group(group)
    if not hexagonal_compatible(ops):raise ValueError("group needs triangular reciprocal lattice")
    if domain=="hexagon":
        k,w,area=hexagonal_bz_quadrature(order)
        sampler="D6 hexagon Gauss"
    elif domain=="primitive":
        k,w,area=primitive_cell_quadrature(order,method=method,shift=shift)
        sampler="primitive "+method
    else:raise ValueError("domain must be hexagon or primitive")
    even=[];odd=[];gaps=[];dist=[];occ=[];active=[]
    for point in k:
        H,V,S=hexagonal_symmetrized_model(point,ops)
        d=local_transition_diagnostic(H,V,S,mu=mu,temperature=temperature,z=z)
        even.append(d["even"]);odd.append(d["odd"])
        gaps.append(d["min_intercluster_gap"])
        dist.append(d["closest_active_resonance_distance"])
        occ.append(d["max_pair_occupation_difference"])
        active.append(d["active_transition_count"])
    e=field_diagnostics(even,w,component=component,top_fraction=top_fraction)
    o=field_diagnostics(odd,w,component=e["component"],top_fraction=top_fraction)
    a,b=e["component"]
    top=e["top_indices"][:hotspot_count]
    hotspots=[dict(k_cart=[float(x) for x in k[i]],
                   k_fractional=[float(x) for x in (INV_BASIS@k[i])],
                   weighted_abs_even=float(w[i]*abs(even[i][a,b])),
                   local_abs_even=float(abs(even[i][a,b])),
                   min_gap=float(gaps[i]),active_resonance_distance=dist[i],
                   max_occupancy_difference=float(occ[i]))
              for i in top]
    return dict(group=group,domain=domain,method=sampler,order=int(order),
                nodes=len(w),bz_area=area,mu=mu,temperature=temperature,z=z,
                component=(a,b),even=e,odd=o,
                min_gap=float(min(gaps)),
                min_active_resonance_distance=float(min(d for d in dist if d is not None))
                    if any(d is not None for d in dist) else None,
                min_occupation_difference_for_active=1e-5,
                max_occupation_difference=float(max(occ)),
                total_active_directed_transition_count=int(sum(active)),
                hotspots=hotspots)
