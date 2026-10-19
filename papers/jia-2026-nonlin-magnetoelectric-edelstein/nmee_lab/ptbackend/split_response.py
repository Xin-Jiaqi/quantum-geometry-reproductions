"""Finite-frequency PT-degenerate EB response with exact Fermi-surface pole isolation.

Conventions follow round27: H_B=Bbar*S, Bbar=-g*muB*B; H_E=eE*x.
Use a fixed orbital basis. A finite-frequency adiabatic pole is NOT automatically
an impurity scattering lifetime; the response coefficient follows the code's
magnetic sign convention. The published Eq.(2) has a sign discrepancy relative
to this convention under current comparison, explicitly retained below.
"""
import numpy as np
from .cumnas import Model, MU_B,IX,IY,IZ,PA,PB
from .liouville_doublet import prepare,resolvent,cm

def split(q,theta,mu,p=Model(),eta_e=1e-4,eta_b=None,spin=IY,obs=IX,beta=1):
    if eta_b is None: eta_b=1.31*eta_e
    Oe=1j*eta_e; Ob=1j*eta_b; Os=Oe+Ob
    d=prepare(q,theta,mu,p,spin,beta)
    Pm,Pp=d['P']; V=d['V'];S=spin
    em,ep=d['em'],d['ep']
    def fd(e):
        if not p.temperature_K: return np.zeros_like(e)
        from scipy.special import expit
        from .cumnas import KB
        f=expit((mu-e)/(KB*p.temperature_K))
        return -f*(1-f)/(KB*p.temperature_K)
    fd_m,fd_p=fd(em),fd(ep)
    vintra=(fd_m[:,None,None]*(Pm@V@Pm)+fd_p[:,None,None]*(Pp@V@Pp))
    vinter=d['drho']-vintra
    # This splitting is an algebraic decomposition, NOT automatically identical
    # to the paper's intrinsic/extrinsic decomposition away from the pole.
    rhoB=resolvent(d,Ob,cm(S,d['rho0']))
    dRhoB0=resolvent(d,Ob,cm(V,rhoB))
    def subtotal(dRho):
        rhoE=1j*resolvent(d,Oe,dRho)
        dRhoB=resolvent(d,Ob,cm(S,dRho))
        rhoEB=resolvent(d,Os,1j*dRhoB+cm(S,rhoE))
        return -p.g*MU_B/p.a_nm**2*np.einsum('ij,nji->n',obs,rhoEB)
    def geom():
        dRhoB=dRhoB0
        rhoEB=resolvent(d,Os,1j*dRhoB)
        return -p.g*MU_B/p.a_nm**2*np.einsum('ij,nji->n',obs,rhoEB)
    parts={'fermi_derivative':subtotal(vintra),'intergroup':subtotal(vinter),'geometry_of_magnetic_response':geom()}
    parts['full']=sum(parts.values())
    return parts

def residue_from_projectors(q,theta,mu,p=Model(),spin=IY,obs=IX,beta=1):
    """Exact 1/eta_E coefficient at fixed nonzero doublet separation.

    gmuB and a^2 as in cumnas. Avoid direct division at a Dirac node.
    Derived from i/[i eta_E] * [spin, f'_a v_a P_a] and the off-diagonal
    static resolvent. No intra-doublet denominators appear.
    """
    data=prepare(q,theta,mu,p,spin,beta)
    pm,pp=data['P']; v=data['V']; e0,e1=data['em'],data['ep']
    from .cumnas import fermi_prime
    if p.temperature_K==0:return np.zeros_like(e0)
    f0=fermi_prime(e0,mu,p.temperature_K)
    f1=fermi_prime(e1,mu,p.temperature_K)
    v0=np.einsum('nij,nji->n',pm,v).real/2
    v1=np.einsum('nij,nji->n',pp,v).real/2
    s01=np.einsum('nij,jk,nkl,li->n',pm,obs,pp,spin).real
    s10=np.einsum('nij,jk,nkl,li->n',pp,obs,pm,spin).real
    assert np.allclose(s01,s10,atol=1e-12)
    return -p.g*MU_B/p.a_nm**2*(2*s01*v0*f0/(e0-e1)+2*s10*v1*f1/(e1-e0))
