"""Round25: finite-T equilibrium Zeeman Ward in a k-dependent Bloch basis.

Fixed orthonormal local orbitals may have a k-dependent unitary Bloch gauge
W(k). Define the anti-Hermitian connection A_a=W^dagger d_a W and
D_a X = d_a X + [A_a,X]. For any smooth k-dependent effective spin vertex
S_b(k), a spatially uniform Zeeman field H_B=H-kappa B_b S_b gives

 chi_JB^{ab}(k,0) = -kappa Tr[(D_a H) Df_H[S_b]]
                    -kappa Tr[f(H) D_a S_b]
                  = -kappa d_a Tr[f(H) S_b].

The second term is the explicit B-dependence of the covariant current
D_a H_B. It is REQUIRED if D_a S_b != 0. This is equilibrium and STATIC;
it neither models finite-q EM contacts nor dynamical dissipative KME.

The synthetic k-dependent effective-spin vertex is not a claim about an
actual PAW/Wannier matrix. No orbital-center embedding or SI calibration.
"""
from __future__ import annotations
import numpy as np
from .causal_response import _fermi
from .core import PAULI
from .degenerate_static_ward import (fermi_divided_difference,
                                      kramers_periodic_model)
from .hexagonal_bz import INV_BASIS


def _validate(H,dH,S,dS,connection,mu,temperature,kappa):
    H=np.asarray(H,complex);dH=np.asarray(dH,complex);S=np.asarray(S,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or len(H)<2:
        raise ValueError("H must be square, N>=2")
    n=len(H)
    if dH.shape!=(2,n,n) or S.shape!=(3,n,n):
        raise ValueError("dH must be (2,N,N), S must be (3,N,N)")
    if dS is None:dS=np.zeros((2,3,n,n),complex)
    else:dS=np.asarray(dS,complex)
    if connection is None:connection=np.zeros((2,n,n),complex)
    else:connection=np.asarray(connection,complex)
    if dS.shape!=(2,3,n,n) or connection.shape!=(2,n,n):
        raise ValueError("dS=(2,3,N,N) and connection=(2,N,N) required")
    if not all(np.isfinite(q) for q in (mu,temperature,kappa)) or temperature<=0:
        raise ValueError("finite mu, kappa and T>0 required")
    for x in (H,*dH,*S,*dS.reshape(-1,n,n)):
        if not np.all(np.isfinite(x)) or not np.allclose(x,x.conj().T,atol=1e-10,rtol=0):
            raise ValueError("H, dH, S, dS must be finite Hermitian")
    for a in connection:
        if not np.all(np.isfinite(a)) or not np.allclose(a,-a.conj().T,atol=1e-10,rtol=0):
            raise ValueError("connection must be finite anti-Hermitian")
    return H,dH,S,dS,connection


def covariant_velocities(H,dH,S,dS,connection):
    """Covariant dH and dS, equivalent to reference physical-frame operators."""
    J=dH+np.stack([a@H-H@a for a in connection])
    DS=dS+np.stack([[a@s-s@a for s in S] for a in connection])
    return J,DS


def fermi_frechet(H,Q,*,mu=.15,temperature=.22):
    """Exact finite-T Frechet action D f_H[Q] including Kramers degeneracy."""
    e,U=np.linalg.eigh(np.asarray(H,complex))
    dd=fermi_divided_difference(e,mu=mu,temperature=temperature)
    Qband=U.conj().T@Q@U
    return U@(dd*Qband)@U.conj().T


def covariant_static_ward(H,dH,S,dS=None,connection=None,*,mu=.15,
                          temperature=.22,kappa=.7):
    """Complete physical static Zeeman current susceptibility and its two terms."""
    H,dH,S,dS,A=_validate(H,dH,S,dS,connection,mu,temperature,kappa)
    J,DS=covariant_velocities(H,dH,S,dS,A)
    e,U=np.linalg.eigh(H);rho=(U*_fermi(e,mu,temperature))@U.conj().T
    frechet=np.stack([fermi_frechet(H,s,mu=mu,temperature=temperature) for s in S])
    population=-kappa*np.einsum("aij,bji->ab",J,frechet,optimize=True)
    contact=-kappa*np.einsum("ij,abji->ab",rho,DS,optimize=True)
    return dict(full_static=population+contact,occupation=population,
                spin_vertex_contact=contact,
                covariant_velocity=J,covariant_spin_gradient=DS,
                spin_density=np.array([np.trace(rho@s) for s in S]),
                status="complete equilibrium static covariant-spin Ward; no dynamic KME")


def finite_zeeman_current_derivative(H,dH,S,dS=None,connection=None,*,mu=.15,
                                    temperature=.22,kappa=.7,step=1e-5):
    """Independent B finite difference including derivative of current operator."""
    H,dH,S,dS,A=_validate(H,dH,S,dS,connection,mu,temperature,kappa)
    if step<=0 or not np.isfinite(step):raise ValueError("positive field step required")
    J,DS=covariant_velocities(H,dH,S,dS,A)
    out=np.zeros((2,3),complex)
    for b in range(3):
        def current(B):
            hb=H-kappa*B*S[b]
            e,U=np.linalg.eigh(hb)
            rho=(U*_fermi(e,mu,temperature))@U.conj().T
            return np.array([np.trace(rho@(J[a]-kappa*B*DS[a,b])) for a in range(2)])
        out[:,b]=(current(step)-current(-step))/(2*step)
    return out


def periodic_spin_embedding(k,*,effective_spin_variation=.16,rotate=True):
    """Synthetic independent periodic H(k), S(k), dH(k), dS(k) and Bloch A(k).

    effective_spin_variation=0 is constant *physical* Pauli spin in reference
    gauge; >0 represents an intentionally k-dependent *effective projected*
    spin vertex. Bloch W(k)=exp(i theta(k) G) rotates all operators.
    """
    k=np.asarray(k,float)
    if k.shape!=(2,) or not np.all(np.isfinite(k)):
        raise ValueError("finite Cartesian k pair required")
    if not np.isfinite(effective_spin_variation):
        raise ValueError("finite spin variation required")
    H,V,S=kramers_periodic_model(k)
    q=INV_BASIS@k
    amp=effective_spin_variation
    g=amp*(np.sin(q[0])+.4*np.cos(q[1]))
    dgq=amp*np.array([np.cos(q[0]),-.4*np.sin(q[1])])
    dg=INV_BASIS.T@dgq
    C=np.stack([np.kron(PAULI[0],s) for s in PAULI])
    S=S+g*C
    dS=np.stack([dg[a]*C for a in range(2)])
    n=len(H)
    if not rotate:
        return H,V,S,dS,np.zeros((2,n,n),complex)
    G=np.kron(PAULI[0],PAULI[1])
    theta=.21*np.sin(q[0])+.14*np.cos(q[1])
    dtheta=INV_BASIS.T@np.array([.21*np.cos(q[0]),-.14*np.sin(q[1])])
    W=np.cos(theta)*np.eye(n)+1j*np.sin(theta)*G
    A=np.stack([1j*dtheta[a]*G for a in range(2)])
    rotate_op=lambda x:W.conj().T@x@W
    h=rotate_op(H);s=np.stack([rotate_op(x) for x in S])
    # Ordinary matrix derivatives in the rotating gauge:
    # d (W^dag X W) = W^dag(dX)W + [X_tilde,A].
    v=np.stack([rotate_op(V[a])+h@A[a]-A[a]@h for a in range(2)])
    ds=np.stack([[rotate_op(dS[a,b])+s[b]@A[a]-A[a]@s[b]
                  for b in range(3)] for a in range(2)])
    return h,v,s,ds,A


def finite_k_spin_density_gradient(k,*,effective_spin_variation=.16,rotate=True,
                                   mu=.15,temperature=.22,kappa=.7,step=1e-5):
    """Independent k finite difference of equilibrium Tr[f(H(k)) S(k)]."""
    k=np.asarray(k,float)
    if k.shape!=(2,) or step<=0 or not np.isfinite(step):
        raise ValueError("finite Cartesian k and positive step required")
    out=np.zeros((2,3),complex)
    for a in range(2):
        dk=np.eye(2)[a]*step
        vals=[]
        for kn in (k+dk,k-dk):
            H,_,S,_,_=periodic_spin_embedding(kn,
                     effective_spin_variation=effective_spin_variation,rotate=rotate)
            e,U=np.linalg.eigh(H)
            rho=(U*_fermi(e,mu,temperature))@U.conj().T
            vals.append(np.array([np.trace(rho@s) for s in S]))
        out[a]=-kappa*(vals[0]-vals[1])/(2*step)
    return out
