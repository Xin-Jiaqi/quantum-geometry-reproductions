"""Stationary electric-response density matrix in explicitly specified CuMnAs RTA.

In a *fixed orbital basis*, H(B)=H0+Bbar*S, rho_eq(B)=f(H(B)).
For fixed static B, solve (i ad_H + C_B) rho_E = D_k rho_eq,
with C_B acting at two rates in the spectral pair-blocks of H(B):
eta_pop for matrix elements in the same PT-descended pair (0,1 or 2,3),
eta_coh for elements connecting different pairs.

Units: eta_pop / eta_coh are energies (hbar / tau), not collision self-energies.
RTA specifies a particular steady-state response; it is NOT the paper's intrinsic
prescription and not a fully microscopic conserving disorder calculation.
The reference equilibrium is allowed to shift with B; this differs from the
fixed-zero-field-occupation adiabatic protocol used in previous rounds.
"""
from dataclasses import replace
import numpy as np
from scipy.special import expit
from .cumnas import Model,hamiltonian,velocity,PA,PB,IX,IZ,MU_B,KB


def divdiff(e,mu,T):
    if T <= 0: raise ValueError('Use finite T: T=0 distribution derivatives require FS integration')
    f=expit((mu-e)/(KB*T))
    fp=-f*(1-f)/(KB*T)
    de=e[:,:,None]-e[:,None,:]
    dif=f[:,:,None]-f[:,None,:]
    out=np.zeros_like(de)
    # T-dependent, numerically stable divided difference around same energies.
    safe=np.abs(de)>max(1e-9,1e-7*KB*T)
    np.divide(dif,de,out=out,where=safe)
    emid=(e[:,:,None]+e[:,None,:])/2
    fm=expit((mu-emid)/(KB*T))
    out=np.where(safe,out,-fm*(1-fm)/(KB*T))
    return out


def steady_chi_from_matrices(H,V,mu,T,obs,a_nm=1.,eta_pop=.01316,eta_coh=.01316,return_parts=False):
    """Same PT-descended pair-block relaxation in a fixed orbital basis.
    Basis-covariant under a k-independent common unitary conjugation of H,V,O.
    """
    e,U=np.linalg.eigh(H)
    udag=np.swapaxes(U.conj(),1,2)
    ve=udag@V@U; oe=udag@obs@U
    ff=divdiff(e,mu,T)
    de=e[:,:,None]-e[:,None,:]
    pair=np.array([0,0,1,1]); is_intra=pair[:,None]==pair[None,:]
    eta=np.where(is_intra,eta_pop,eta_coh)
    if eta_pop<=0 or eta_coh<=0:raise ValueError("rates must be positive")
    rhoe=ff*ve/(eta[None,:,:]+1j*de)
    val=np.einsum("nij,nji->n",oe,rhoe).real/(a_nm**2)
    if return_parts:
        p1=np.einsum("nij,nji->n",oe,rhoe*is_intra).real/(a_nm**2)
        p2=np.einsum("nij,nji->n",oe,rhoe*(~is_intra)).real/(a_nm**2)
        return val,p1,p2
    return val

def steady_chi(q,theta,mu,par=Model(),Bbar=0.,spin=IX,obs=PA@IZ,beta=0,
               eta_pop=.01316,eta_coh=.01316, return_parts=False):
    """Compute chi_E(Bbar) per electric field in specified RTA.
    Normalization omits -g*muB; an EB response is obtained by derivative w.r.t Bbar.
    Blocks continue the two zero-B PT multiplets; requires nonclosing interpair gap.
    """
    q=np.atleast_2d(q)
    H=hamiltonian(q,theta,par)+Bbar*spin
    V=velocity(q,theta,beta,par)
    e,U=np.linalg.eigh(H)
    gap=np.min(e[:,2]-e[:,1])
    if gap<=5*abs(Bbar):raise ValueError(f'PT pair mixing / near gap closing, gap={gap}')
    return steady_chi_from_matrices(H,V,mu,par.temperature_K,obs,par.a_nm,eta_pop,eta_coh,return_parts)

def nm_response(q,theta,mu,par=Model(),eta_pop=.01316,eta_coh=.01316,
                step_bar=1e-6,spin=IX,obs=PA@IZ,beta=0,return_parts=False):
    """d chi_E/d Bbar *(-g mu_B), per physical B in T.
    Default spin=IX and obs=PA@IZ represent Fig.1c.
    """
    d0=steady_chi(q,theta,mu,par,step_bar,spin,obs,beta,eta_pop,eta_coh,return_parts)
    d1=steady_chi(q,theta,mu,par,-step_bar,spin,obs,beta,eta_pop,eta_coh,return_parts)
    if return_parts:
        return np.stack([(-par.g*MU_B)*(x-y)/(2*step_bar) for x,y in zip(d0,d1)])
    return -par.g*MU_B*(d0-d1)/(2*step_bar)


def integrate(q,*args,**kwargs):
    return nm_response(q,*args,**kwargs).mean(axis=-1)
