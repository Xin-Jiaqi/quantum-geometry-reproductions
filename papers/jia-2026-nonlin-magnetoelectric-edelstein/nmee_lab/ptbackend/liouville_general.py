"""Finite-frequency EB density-matrix response for arbitrary Hermitian H(k) with degeneracies.

FIXED ORBITAL GAUGE ONLY: V=∂kH, S(k) constant. For a k-dependent
orbital gauge, replace every derivative by the *orbital connection*
covariant derivative; supplied velocity and D_k spin alone are not enough
unless the gauge-covariant form of the derivative is consistently used.

Retarded energy variables Ω=ℏω+iη in eV; ωE and ωB are drive frequencies.
Returns response in μB magnetic moment units using convention H_B = -gμB B S,
H_E = e E x, with eE in eV/nm and spin observable dimensionless.
"""
import numpy as np
from .cumnas import MU_B, KB
from scipy.special import expit

def comm(a,b):return a@b-b@a

def spec_data(H,V,mu,T):
    e,U=np.linalg.eigh(H)
    Ud=np.swapaxes(U.conj(),-1,-2)
    Ve=Ud@V@U
    if T>0:
        f=expit((mu-e)/(KB*T)); fp= -f*(1-f)/(KB*T)
    else:
        f=(e<mu).astype(float);fp=np.zeros_like(e)
    fdiv=np.empty((*e.shape[:-1],len(e[0]),len(e[0])))
    de=e[..., :,None]-e[...,None,:]
    num=f[..., :,None]-f[...,None,:]
    # near degeneracy, use analytical f' at average energy to avoid catastrophic cancellation
    close=np.abs(de)<1.e-7
    np.divide(num,de,out=fdiv,where=~close)
    if T>0:
        mid=(e[..., :,None]+e[..., None,:])/2
        fm=expit((mu-mid)/(KB*T));
        np.copyto(fdiv,-fm*(1-fm)/(KB*T),where=close)
    else:
        np.copyto(fdiv,fp[..., :,None],where=close)
    drho=U@(fdiv*Ve)@Ud
    rho0=U@(f[..., :,None]*np.eye(len(e[0])))@Ud
    return e,U,Ud,drho,rho0

def R(e,U,Ud,X,Omega):
    Xe=Ud@X@U
    delta=e[..., :,None]-e[...,None,:]
    return U@(Xe/(Omega-delta))@Ud

def eb_single(H,V,S,O,mu,T,g,OmegaE,OmegaB, a_nm=1.,dSpin=None, return_parts=False):
    H=np.asarray(H,dtype=complex);V=np.asarray(V,dtype=complex)
    if H.ndim==2:H=H[None,:,:];V=V[None,:,:]
    e,U,Ud,drho,rho0=spec_data(H,V,mu,T)
    rhoE=1j*R(e,U,Ud,drho,OmegaE)
    rhoB=R(e,U,Ud,comm(S,rho0),OmegaB)
    drhoB=R(e,U,Ud,comm(V,rhoB)+comm(S,drho)+(comm(dSpin,rho0) if dSpin is not None else 0),OmegaB)
    rhoEB=R(e,U,Ud,1j*drhoB+comm(S,rhoE),OmegaE+OmegaB)
    val=(-g*MU_B/a_nm**2)*np.einsum('ij,nji->n',O,rhoEB)
    return (val,rhoEB) if return_parts else val

def eb_cos(H,V,S,O,mu,T,g,omegaE,omegaB,eta,a_nm=1.):
    ret=0.
    for se in (-1,1):
        for sb in (-1,1):
            ret+=eb_single(H,V,S,O,mu,T,g,se*omegaE+1j*eta,sb*omegaB+1j*eta,a_nm)/4
    return ret
