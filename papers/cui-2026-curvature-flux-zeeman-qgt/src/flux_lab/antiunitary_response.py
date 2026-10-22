"""Microscopic 2D antiunitary symmetry oracles for Zeeman interband response.

Physical spin is Pauli, not an orbital pseudospin. Toy models are synthetic.
A=U K, k'=Rk, U sigma_b* U^dag=eps_b sigma_b.
T_ab(Rk)=-R_ac eps_b T_cb(k)*; nondegenerate band matching.
chi_JB(Rk,z)=R chi_JB(k,-z) diag(eps)
alpha_ME(Rk,z)=-diag(eps) alpha_ME(k,-z) R.T.
For z=omega+i eta, -z is advanced; -omega+i eta is NOT -z.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .core import PAULI
from .causal_response import _data, interband_kubo, interband_liouville

@dataclass(frozen=True)
class Antiunitary:
    name:str
    R:np.ndarray
    U:np.ndarray
    spin_parity:np.ndarray

def examples():
    return {
        "TR":Antiunitary("1'",-np.eye(2),1j*PAULI[1],np.array([-1.,-1.,-1.])),
        "C2zT":Antiunitary("2'",np.eye(2),PAULI[0],np.array([1.,1.,-1.])),
        "mxT":Antiunitary("m'",np.diag([1.,-1.]),PAULI[2],np.array([-1.,1.,1.]))
    }

def toy_bloch(k,kind="TR"):
    x,y=np.asarray(k,float)
    if kind=="TR":
        d=np.array([.81*np.sin(x),.69*np.sin(y),.27*np.sin(x+y)])
        derivative=np.array([[.81*np.cos(x),0,.27*np.cos(x+y)],
                             [0,.69*np.cos(y),.27*np.cos(x+y)]])
        scalar=.12*np.cos(x)+.07*np.cos(y)
        dscalar=np.array([-.12*np.sin(x),-.07*np.sin(y)])
    elif kind=="C2zT":
        d=np.array([.71+.31*np.sin(x)+.12*np.cos(y),
                    .35*np.sin(y)+.16*np.cos(x),0.])
        derivative=np.array([[.31*np.cos(x),-.16*np.sin(x),0.],
                             [-.12*np.sin(y),.35*np.cos(y),0.]])
        scalar=.13*np.cos(x)+.04*np.sin(y)
        dscalar=np.array([-.13*np.sin(x),.04*np.cos(y)])
    elif kind=="mxT":
        d=np.array([.8*np.sin(y),.34+.3*np.cos(x),.59+.13*np.cos(y)])
        derivative=np.array([[0.,-.3*np.sin(x),0.],
                             [.8*np.cos(y),0.,-.13*np.sin(y)]])
        scalar=.11*np.cos(x)+.02*np.cos(y)
        dscalar=np.array([-.11*np.sin(x),-.02*np.sin(y)])
    else:raise ValueError("kind: TR, C2zT, or mxT")
    H=np.einsum('a,aij->ij',d,PAULI)+scalar*np.eye(2)
    V=np.einsum('ba,aij->bij',derivative,PAULI)+dscalar[:,None,None]*np.eye(2)
    return H,V,PAULI.copy()

def _checks(spec):
    if spec.R.shape!=(2,2) or spec.U.shape!=(2,2) or spec.spin_parity.shape!=(3,):
        raise ValueError("invalid antiunitary spinor mapping")
    if not np.allclose(spec.R@spec.R.T,np.eye(2),atol=1e-12):
        raise ValueError("momentum mapping must be orthogonal")
    if not np.allclose(spec.U@spec.U.conj().T,np.eye(2),atol=1e-12):
        raise ValueError("antiunitary spinor U must be unitary")
    mapped=np.stack([spec.U@x.conj()@spec.U.conj().T for x in PAULI])
    if not np.allclose(mapped,spec.spin_parity[:,None,None]*PAULI,atol=1e-12):
        raise ValueError("spin parity incompatible with U")

def operator_residual(kind,k):
    spec=examples()[kind];_checks(spec)
    H,V,S=toy_bloch(k,kind);H1,V1,_=toy_bloch(spec.R@k,kind)
    trans=lambda A:spec.U@A.conj()@spec.U.conj().T
    vp=np.einsum('ab,bij->aij',spec.R,np.stack([trans(q) for q in V]))
    return dict(H=float(np.max(abs(H1-trans(H)))),
                velocity=float(np.max(abs(V1-vp))),
                spin=float(np.max(abs(np.stack([trans(q) for q in S])-
                                            spec.spin_parity[:,None,None]*S))))

def qgt_band_pair(h,v,s,n=0,m=1):
    D=_data(h,v,s)
    if n==m:raise ValueError('interband n != m')
    return np.outer(D['r'][:,n,m],D['s'][:,m,n])

def qgt_symmetry_residual(kind,k,n=0,m=1):
    spec=examples()[kind]
    h,v,s=toy_bloch(k,kind);hh,vv,ss=toy_bloch(spec.R@k,kind)
    A=qgt_band_pair(h,v,s,n,m);B=qgt_band_pair(hh,vv,ss,n,m)
    prediction=-np.einsum('ac,cb,b->ab',spec.R,A.conj(),spec.spin_parity)
    return float(np.max(abs(B-prediction)))

def kubo_symmetry_residual(kind,k,z=.07+.014j,mu=.1,temperature=.16,method='kubo'):
    spec=examples()[kind];h,v,s=toy_bloch(k,kind)
    hp,vp,sp=toy_bloch(spec.R@k,kind)
    opts=dict(z=z,mu=mu,temperature=temperature)
    fn=interband_kubo if method=='kubo' else interband_liouville if method=='liouville' else None
    if fn is None:raise ValueError('method must be kubo or liouville')
    chi1,alpha1=fn(hp,vp,sp,**opts)
    chi0,alpha0=fn(h,v,s,**{**opts,'z':-z})
    predchi=spec.R@chi0@np.diag(spec.spin_parity)
    predalpha=-np.diag(spec.spin_parity)@alpha0@spec.R.T
    return dict(chi=float(np.max(abs(chi1-predchi))),
                alpha=float(np.max(abs(alpha1-predalpha))),
                chi_at_partner=chi1,alpha_at_partner=alpha1,
                predicted_chi=predchi,predicted_alpha=predalpha)

def pair_frequency_parities(kind,k,z=.027,mu=.1,temperature=.16):
    if complex(z).imag!=0:raise ValueError('clean real z required')
    spec=examples()[kind]
    h,v,s=toy_bloch(k,kind);hp,vp,sp=toy_bloch(spec.R@k,kind)
    opts=dict(z=float(z),mu=mu,temperature=temperature)
    pos=interband_kubo(h,v,s,**opts)[0]+interband_kubo(hp,vp,sp,**opts)[0]
    neg=interband_kubo(h,v,s,**{**opts,'z':-z})[0]+interband_kubo(hp,vp,sp,**{**opts,'z':-z})[0]
    return dict(positive=pos,negative=neg,odd=(pos-neg)/2,even=(pos+neg)/2,
                parity_matrix=np.diag(spec.R)[:,None]*spec.spin_parity[None,:])

def wrong_retarded_branch_error(kind,k,omega=.08,eta=.025):
    spec=examples()[kind]
    h,v,s=toy_bloch(k,kind);hp,vp,sp=toy_bloch(spec.R@k,kind)
    z=complex(omega,eta);opts=dict(mu=.1,temperature=.16)
    target=interband_kubo(hp,vp,sp,z=z,**opts)[0]
    wrong=spec.R@interband_kubo(h,v,s,z=complex(-omega,eta),**opts)[0]@np.diag(spec.spin_parity)
    return float(np.max(abs(target-wrong)))

def toy_spin_orbital(k,kind="TR",orbital_split=1.35,mixing=.21):
    """Four-band real physical spin and coupled orbital tau, not a material."""
    H,V,S=toy_bloch(k,kind)
    tauz=PAULI[2];taux=PAULI[0];I=np.eye(2)
    H4=np.kron(I,H)+orbital_split*np.kron(tauz,I)+mixing*np.kron(taux,H)
    V4=np.stack([np.kron(I,v)+mixing*np.kron(taux,v) for v in V])
    S4=np.stack([np.kron(I,s) for s in S])
    return H4,V4,S4

def fourband_antiunitary_residual(kind,k,z=.07+.014j):
    spec=examples()[kind];U4=np.kron(np.eye(2),spec.U)
    H,V,S=toy_spin_orbital(k,kind)
    HH,VV,SS=toy_spin_orbital(spec.R@k,kind)
    trans=lambda A:U4@A.conj()@U4.conj().T
    residual_H=np.max(abs(HH-trans(H)))
    residual_V=np.max(abs(VV-np.einsum('ab,bij->aij',spec.R,
                                       np.stack([trans(x) for x in V]))))
    opts=dict(z=z,mu=.1,temperature=.16)
    chi0,alpha0=interband_kubo(H,V,S,**{**opts,'z':-z})
    chi1,alpha1=interband_kubo(HH,VV,SS,**opts)
    return dict(H=float(residual_H),velocity=float(residual_V),
                chi=float(np.max(abs(chi1-spec.R@chi0@np.diag(spec.spin_parity)))),
                alpha=float(np.max(abs(alpha1+np.diag(spec.spin_parity)@alpha0@spec.R.T))))
