"""Periodic triangular-lattice model and symmetry-adapted hexagonal BZ quadrature.

Coordinates are Cartesian k with reciprocal primitive vectors 2*pi*B.
The first BZ is the Wigner-Seitz hexagon of radius 2*pi/sqrt(3).
The quadrature triangulates the hexagon into six equal sectors, then averages
the rule with its mirror image so the weighted node set is exactly D6
symmetric. C3/C6 models are periodic in triangular reciprocal coordinates.
All response outputs are in Cartesian axes.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss
from .core import PAULI
from .group_microscopic import build_group,momentum_map,spin_axial
from .nonabelian_response import cluster_interband_kubo

I=np.eye(2,dtype=complex)
BASIS=np.array([[1.,.5],[0.,np.sqrt(3.)/2]])
INV_BASIS=np.linalg.inv(BASIS)
SIGMA=np.stack([np.kron(I,s) for s in PAULI])


def fractional_action(g,*,integer_tolerance=1e-8):
    """Represent Cartesian antiunitary momentum map in triangular reciprocal coordinates."""
    raw=INV_BASIS@momentum_map(g)@BASIS
    A=np.rint(raw)
    if not np.allclose(raw,A,rtol=0,atol=integer_tolerance):
        raise ValueError("operation does not preserve the triangular reciprocal lattice")
    if abs(round(np.linalg.det(A)))!=1:
        raise ValueError("fractional reciprocal action must be unimodular")
    return A.astype(int)


def hexagonal_compatible(group):
    try:
        ops=build_group(group) if isinstance(group,str) else group
        for g in ops:fractional_action(g)
        return True
    except ValueError:
        return False


def _seed(q):
    """Generic periodic four-band orbital x physical-spin seed and analytic gradient."""
    q=np.asarray(q,float)
    if q.shape!=(2,) or not np.all(np.isfinite(q)):raise ValueError("q must be finite length 2")
    tau=np.stack([np.kron(t,I) for t in PAULI])
    H=2.35*tau[2]+.23*tau[0]
    V=np.zeros((2,4,4),complex)
    terms=[
      ((1,0),"cos",.37*tau[2]),
      ((0,1),"cos",-.29*tau[2]),
      ((1,1),"cos",.21*tau[0]),
      ((1,-1),"sin",.31*np.kron(PAULI[0],PAULI[0])),
      ((1,0),"sin",.26*np.kron(PAULI[0],PAULI[1])),
      ((0,1),"sin",-.19*np.kron(PAULI[0],PAULI[2])),
      ((1,1),"sin",.17*np.kron(PAULI[1],PAULI[0])),
      ((2,1),"cos",.13*np.kron(PAULI[1],PAULI[2])),
      ((1,2),"sin",.11*np.kron(PAULI[2],PAULI[1])),
      ((2,-1),"cos",.09*np.kron(PAULI[2],PAULI[0])),
    ]
    for n,kind,A in terms:
        n=np.asarray(n,float);phase=float(n@q)
        if kind=="sin":f,df=np.sin(phase),np.cos(phase)
        else:f,df=np.cos(phase),-np.sin(phase)
        H=H+f*A
        V+=n[:,None,None]*df*A
    return H,V


def hexagonal_symmetrized_model(k_cart,group):
    """Return H(k), Cartesian dH/dk, physical spin for a hexagonal-lattice MPG."""
    k=np.asarray(k_cart,float)
    if k.shape!=(2,) or not np.all(np.isfinite(k)):raise ValueError("k_cart must be finite length 2")
    ops=build_group(group) if isinstance(group,str) else group
    actions=[fractional_action(g).astype(float) for g in ops]
    q=INV_BASIS@k
    H=np.zeros((4,4),complex);Vq=np.zeros((2,4,4),complex)
    for g,A in zip(ops,actions):
        q0=np.linalg.solve(A,q)
        h,v=_seed(q0)
        U4=np.kron(I,g.U)
        trans=lambda X:U4@(X.conj() if g.anti else X)@U4.conj().T
        H+=trans(h)
        Vq+=np.einsum("ab,bij->aij",np.linalg.inv(A).T,
                      np.stack([trans(x) for x in v]))
    H/=len(ops);Vq/=len(ops)
    V=np.einsum("ab,bij->aij",INV_BASIS.T,Vq)
    return H,V,SIGMA.copy()


def hexagonal_bz_quadrature(order=6):
    """Symmetry-adapted Gauss rule over the first hexagonal Wigner-Seitz BZ.

    Returns k_cart (12*order**2,2), normalized area weights summing to one,
    and the exact dimensionless BZ area 2*sqrt(3)*pi**2.
    """
    if not isinstance(order,(int,np.integer)) or order<2:
        raise ValueError("order must be integer >=2")
    x,w=leggauss(order);u=(x+1)/2;wu=w/2
    radius=2*np.pi/np.sqrt(3.)
    angles=np.pi/6+np.arange(6)*np.pi/3
    vertices=radius*np.stack([np.cos(angles),np.sin(angles)],axis=1)
    points=[];weights=[]
    for i in range(6):
        a=vertices[i];b=vertices[(i+1)%6]
        jac=abs(np.linalg.det(np.column_stack([a,b])))
        for iu,uu in enumerate(u):
            for iv,vv in enumerate(u):
                points.append(uu*a+(1-uu)*vv*b)
                weights.append(jac*(1-uu)*wu[iu]*wu[iv])
    points=np.asarray(points);weights=np.asarray(weights)
    area=2*np.sqrt(3.)*np.pi**2
    if abs(np.sum(weights)-area)>2e-12*area:
        raise AssertionError("hexagonal quadrature weights do not cover BZ area")
    # The triangle map is C6 symmetric but not individually mirror symmetric.
    # Average with one reflection; C6 then generates all six mirror operations.
    mirror=np.diag([-1.,1.])
    points=np.concatenate([points,points@mirror.T],axis=0)
    weights=np.concatenate([weights/2,weights/2])
    return points,weights/area,float(area)


def hexagonal_local_residual(group,k_cart=np.array([.47,.31]),z=.073+.006j,
                             mu=.13,temperature=.17):
    """Pointwise H/velocity and full complex-frequency Kubo covariance check."""
    ops=build_group(group)
    if not hexagonal_compatible(ops):
        raise ValueError("group is incompatible with triangular reciprocal lattice")
    k=np.asarray(k_cart,float);H,V,S=hexagonal_symmetrized_model(k,ops)
    errors=np.zeros(4)
    for g in ops:
        Q=momentum_map(g);L=spin_axial(g)
        Hp,Vp,Sp=hexagonal_symmetrized_model(Q@k,ops)
        U4=np.kron(I,g.U)
        trans=lambda X:U4@(X.conj() if g.anti else X)@U4.conj().T
        errors[0]=max(errors[0],np.max(abs(Hp-trans(H))))
        predv=np.einsum("ab,bij->aij",Q,np.stack([trans(x) for x in V]))
        errors[1]=max(errors[1],np.max(abs(Vp-predv)))
        chi0,alpha0=cluster_interband_kubo(H,V,S,z=-z if g.anti else z,mu=mu,temperature=temperature)
        chi1,alpha1=cluster_interband_kubo(Hp,Vp,Sp,z=z,mu=mu,temperature=temperature)
        errors[2]=max(errors[2],np.max(abs(chi1-Q@chi0@L.T)))
        errors[3]=max(errors[3],np.max(abs(alpha1-(-1 if g.anti else 1)*L@alpha0@Q.T)))
    return dict(H=float(errors[0]),velocity=float(errors[1]),
                chi=float(errors[2]),alpha=float(errors[3]),order=len(ops))


def hexagonal_bz_response(group,*,order=5,z=.071+.004j,mu=.13,temperature=.17):
    """Integrate even/odd causal response over the true first hexagonal BZ.

    The negative-frequency branch is -z (advanced when Im z>0), not
    -Re(z)+i Im(z). The returned integral is normalized by BZ area.
    """
    ops=build_group(group)
    if not hexagonal_compatible(ops):
        raise ValueError("group is incompatible with triangular reciprocal lattice")
    points,weights,area=hexagonal_bz_quadrature(order)
    plus=np.zeros((2,3),complex);minus=np.zeros((2,3),complex)
    for k,w in zip(points,weights):
        H,V,S=hexagonal_symmetrized_model(k,ops)
        a=cluster_interband_kubo(H,V,S,z=z,mu=mu,temperature=temperature)[0]
        b=cluster_interband_kubo(H,V,S,z=-z,mu=mu,temperature=temperature)[0]
        plus+=w*a;minus+=w*b
    even=(plus+minus)/2;odd=(plus-minus)/2
    return dict(even=even,odd=odd,positive_frequency=plus,negative_frequency=minus,
                points=points,weights=weights,bz_area=area,order=order)
