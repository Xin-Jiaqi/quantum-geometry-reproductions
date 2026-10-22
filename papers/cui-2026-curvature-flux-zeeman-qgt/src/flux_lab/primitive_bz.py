"""Independent primitive-reciprocal-cell integration for the Cui Flux toy model.

A primitive reciprocal parallelogram and the first Wigner-Seitz hexagon are
different fundamental domains of the same triangular reciprocal torus. Exact
integrals of periodic functions must agree; numerical quadrature need not.

We integrate Cartesian physical tensors using k=B q, q in [-pi,pi)^2.
Both the normalized primitive-cell measure and the hexagonal integration
are area-average d2k / Area_BZ. We deliberately DO NOT project the primitive
result into a symmetry-invariant tensor subspace: doing so would hide errors.

The integration routines only support the existing clean *interband cluster*
synthetic model. Neither a material nor metallic intraband conductivity.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss
from .hexagonal_bz import BASIS,hexagonal_symmetrized_model,hexagonal_compatible
from .group_microscopic import build_group
from .nonabelian_response import cluster_interband_kubo


def primitive_cell_quadrature(order, *,method="gauss", shift=(0.,0.)):
    """Quadrature on reciprocal parallelogram: normalized weights sum to 1.

    method='gauss': tensor-product Gauss-Legendre over fractional q axes.
    method='midpoint': shifted periodic reciprocal-coordinate midpoint grid,
    shift in fractional cell units. These are different integration protocols.
    """
    if not isinstance(order,(int,np.integer)) or order<2:
        raise ValueError("order must be an integer >=2")
    if method=="gauss":
        if not np.allclose(shift,(0.,0.),atol=0):
            raise ValueError("shift only applies to periodic midpoint grids")
        u,w=leggauss(int(order))
        nodes=np.pi*u;weights=w/2
    elif method=="midpoint":
        shift=np.asarray(shift,float)
        if shift.shape!=(2,) or not np.all(np.isfinite(shift)):
            raise ValueError("shift must be a finite pair")
        if np.any(np.abs(shift)>=.5):
            raise ValueError("fractional midpoint shift must have |shift|<0.5")
        nodes=np.arange(order,dtype=float)
        # Nodes are [−pi,pi) modulo periodic identification, shifted in
        # index units; independent offsets are applied per coordinate.
        qs=-np.pi+2*np.pi*(nodes+.5)/order
        q1=qs+2*np.pi*shift[0]/order
        q2=qs+2*np.pi*shift[1]/order
        q=np.array([(x,y) for x in q1 for y in q2],float)
        k=q@BASIS.T
        w=np.full(order*order,1./(order*order))
        return k,w,float(4*np.pi**2*abs(np.linalg.det(BASIS)))
    else:
        raise ValueError("method must be 'gauss' or 'midpoint'")
    q=np.array([(x,y) for x in nodes for y in nodes],float)
    weights=np.array([wx*wy for wx in weights for wy in weights],float)
    k=q@BASIS.T
    return k,weights,float(4*np.pi**2*abs(np.linalg.det(BASIS)))


def primitive_bz_response(group, *,order=10,method="gauss",shift=(0.,0.),
                          z=.071+.004j,mu=.13,temperature=.17):
    """Unprojected independent primitive-cell integration.

    Evaluates BOTH complex-frequency signs before even/odd decomposition.
    Group operation and reciprocal normalization match the hex BZ solver;
    sampling domains and quadrature constructions are independent.
    """
    ops=build_group(group) if isinstance(group,str) else group
    if not hexagonal_compatible(ops):
        raise ValueError("group incompatible with triangular reciprocal lattice")
    points,weights,area=primitive_cell_quadrature(order,method=method,shift=shift)
    pos=np.zeros((2,3),complex)
    neg=np.zeros((2,3),complex)
    for k,w in zip(points,weights):
        H,V,S=hexagonal_symmetrized_model(k,ops)
        pos+=w*cluster_interband_kubo(H,V,S,z=z,mu=mu,temperature=temperature)[0]
        neg+=w*cluster_interband_kubo(H,V,S,z=-z,mu=mu,temperature=temperature)[0]
    return {"even":(pos+neg)/2,"odd":(pos-neg)/2,
            "positive_frequency":pos,"negative_frequency":neg,
            "bz_area":area,"order":order,"method":method,"shift":tuple(shift),
            "number_points":len(points),"weight_sum":float(weights.sum())}
