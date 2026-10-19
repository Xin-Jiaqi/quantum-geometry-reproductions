"""Finite-frequency clean EB response in a fixed-orbital CuMnAs basis.
Retains exact twofold PT degeneracy via spectral projectors, not eigenvector labels.
Gauge: fixed orbital basis, no extra k-dependent orbital-position connection.
This is *full finite-frequency* response, not automatically Eq.(9)'s intrinsic piece.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.special import expit
from .cumnas import Model, terms, hamiltonian, velocity, energies, TX,TZ,IX,IY,IZ,PA,PB,KB,MU_B

EYE=np.eye(4,dtype=complex)

def cm(a,b):return a@b-b@a

def prepare(q,theta,mu,p=Model(),spin=IX,beta=0):
    """Prepare vectorized exact projector data, including Frechet derivative of f(H)."""
    q=np.atleast_2d(np.asarray(q,float))
    H=hamiltonian(q,theta,p); V=velocity(q,theta,beta,p)
    em,ep,d=energies(q,theta,p)
    if np.any(d<1e-10): raise ValueError('At exact Dirac touching the two-group projector split is singular. Use finite q-distance or a full 4x4 resolvent.')
    D=H-(em+ep)[:,None,None]*.5*EYE
    Pm=(EYE-D/d[:,None,None])*.5; Pp=(EYE+D/d[:,None,None])*.5
    if p.temperature_K:
        fm=expit((mu-em)/(KB*p.temperature_K)); fp=expit((mu-ep)/(KB*p.temperature_K))
        fpm=-fm*(1-fm)/(KB*p.temperature_K); fpp=-fp*(1-fp)/(KB*p.temperature_K)
    else:
        fm=(em<mu).astype(float); fp=(ep<mu).astype(float);fpm=np.zeros_like(fm);fpp=np.zeros_like(fp)
    rho0=fm[:,None,None]*Pm+fp[:,None,None]*Pp
    cross=(fm-fp)/(em-ep)
    dk_rho0=(fpm[:,None,None]*(Pm@V@Pm)+fpp[:,None,None]*(Pp@V@Pp)+cross[:,None,None]*(Pm@V@Pp+Pp@V@Pm))
    return dict(H=H,V=V,P=(Pm,Pp),em=em,ep=ep,d=d,rho0=rho0,drho=dk_rho0, spin=spin)

def resolvent(data, omega_eV, X):
    """(omega-L_H)^(-1) X with exact twofold spectral projectors.
    Nonzero complex omega ensures all channels well-defined.
    """
    Pm,Pp=data['P'];em,ep=data['em'],data['ep']
    pp=Pm@X@Pm+Pp@X@Pp
    pm=Pm@X@Pp;mp=Pp@X@Pm
    return pp/omega_eV +pm/(omega_eV+(ep-em))[:,None,None]+mp/(omega_eV-(ep-em))[:,None,None]

def eb_kernel(q,theta,mu, p=Model(), eomega=1e-3+1e-4j, bomega=1.31e-3+1e-4j, observation=PA@IZ, spin=IX, beta=0, include_fprime=True):
    """rho_EB per E*Bbar in natural units.
    Output Gamma = (-g mu_B/a^2) Tr(O rho_EB) for response in mu_B units,
    with velocity carrying a_nm; 2D BZ average is separate.
    """
    data=prepare(q,theta,mu,p,spin,beta)
    drho=data['drho'] if include_fprime else None
    if not include_fprime:
        fm=fermi=data['rho0']; Pm,Pp=data['P'];d=data['d'];V=data['V'];
        # Freeze occupation derivatives in electric vertex; not a conserving physical response
        em,ep=data['em'],data['ep']; f1=np.real(np.einsum('nij,nji->n',Pm,fm))/2; f2=np.real(np.einsum('nij,nji->n',Pp,fm))/2
        c=(f1-f2)/(em-ep)
        drho=c[:,None,None]*(Pm@V@Pp+Pp@V@Pm)
    rhoE=1j*resolvent(data,eomega,drho)
    rhoB=resolvent(data,bomega,cm(spin,data['rho0']))
    drhoB=resolvent(data,bomega,cm(data['V'],rhoB)+cm(spin,drho))
    rhoEB=resolvent(data,eomega+bomega,1j*drhoB+cm(spin,rhoE))
    trace=np.einsum('ij,nji->n',observation,rhoEB)
    return (-p.g*MU_B/(p.a_nm**2))*trace, data

def sym_pair_kernel(q,theta,mu,p=Model(),omegaE=1e-3,omegaB=1.31e-3, eta=1e-5,obs=PA@IZ,**kw):
    """Real-field cos drive: sum frequencies ± omega (with +i eta retarded)."""
    ans=0
    for se in (1.,-1.):
        for sb in (1.,-1.):
            res,_=eb_kernel(q,theta,mu,p,complex(se*omegaE,eta),complex(sb*omegaB,eta),obs,**kw)
            ans=ans+res/4.
    return ans
