"""Particle-conserving, spatially uniform BGK for CuMnAs (phenomenological).

Physical protocol: thermalize at fixed weak Zeeman B, then apply homogeneous
DC electric field.  Same two-band-pair dephasing as round29; allow k-dependent
intra-pair population relaxation eta_pop(k).

Dissipator: C[rho](k)= C0[rho](k)-eta_pop(k)*W_B(k)*delta_mu_C[rho],
 W_B = d f(H_B)/d mu = - f'(H_B),
 delta_mu_C = int Tr[C0 rho] / int Tr[eta_pop W_B].
This yields int Tr C=0 for arbitrary rho.  Steady solution of
 i[H_B,rho]+C[rho] = d_k f(H_B) satisfies neutrality int Tr rho=0.

For finite, discrete quadrature exact uniform source charge neutrality may
be imperfect: an optional source-charge centering test reveals the mesh error.

NOT a microscopic impurity self-energy, vertex correction, nor the PRL's
unique intrinsic NMEE; eta_pop anisotropy is illustrative.
"""
import numpy as np
from scipy.special import expit
from .cumnas import Model,hamiltonian,velocity,KB,IX,IY,IZ,PA,MU_B,points
from .steady_relaxation import divdiff


def tau_profile(q, theta, alpha=0.0):
    """Positive dimensionless tau(k)/tau0. Illustration of PT-allowed
    momentum anisotropy; not inferred from microscopic CuMnAs scattering.
    """
    q = np.atleast_2d(q)
    if abs(alpha)>=0.95: raise ValueError('Need positive and bounded tau')
    return 1.0 + alpha*(np.cos(theta)*np.sin(q[:,1])-np.sin(theta)*np.sin(q[:,0]))


def moments(q,theta,mu,Bbar=0.,p=Model(),spin=IX,beta=0,obs=PA@IZ,
            eta_p0=0.013164239138,eta_phi=0.013164239138,
            alpha=0., return_detail=False):
    q=np.atleast_2d(q)
    H=hamiltonian(q,theta,p)+Bbar*spin
    V=velocity(q,theta,beta,p)
    e,U=np.linalg.eigh(H)
    de=e[:,:,None]-e[:,None,:]
    ue=np.swapaxes(U.conj(),1,2)
    Ve=ue@V@U
    Oe=ue@obs@U
    ff=divdiff(e,mu,p.temperature_K)
    tau=tau_profile(q,theta,alpha)
    eta_p=eta_p0/tau
    pair=np.array([0,0,1,1]); intra=(pair[:,None]==pair[None,:])
    enm=np.where(intra[None,:,:],eta_p[:,None,None],eta_phi)
    rho=ff*Ve/(enm+1j*de)
    rta=np.einsum('nij,nji->n',Oe,rho).real/(p.a_nm**2)
    # X0 local occupation source and charge correction moments; no need to
    # construct W matrix, since diagonal in exact H(B) eigenbasis.
    f=expit((mu-e)/(KB*p.temperature_K))
    W=f*(1-f)/(KB*p.temperature_K)
    Vdiag=np.einsum('nii->ni',Ve).real
    Odiag=np.einsum('nii->ni',Oe).real
    source_charge=np.sum((-W)*Vdiag,axis=1)   # Trace d_k f(H)
    X_charge=source_charge/eta_p
    obsW=(Odiag*W).sum(axis=1)/(p.a_nm**2)
    Wtrace=W.sum(axis=1)
    # Uniform fixed-particle-number ensemble constraint ∫ Tr delta_rho=0.
    shift_mu=-X_charge.mean()/Wtrace.mean()
    bgk=rta+shift_mu*obsW
    res=dict(rta=float(rta.mean()),bgk=float(bgk.mean()),correction=float(shift_mu*obsW.mean()),
             delta_mu=float(shift_mu),n_rta=float(X_charge.mean()),
             n_bgk=float(X_charge.mean()+shift_mu*Wtrace.mean()),
             source_total=float(source_charge.mean()),
             observable_W=float(obsW.mean()), W_total=float(Wtrace.mean()),
             # Collision-charge residual is mathematically equal to source
             collision_source_residual=float(source_charge.mean()))
    if return_detail:
        res.update(k_rta=rta,k_bgk=bgk,k_X_charge=X_charge,k_Wtrace=Wtrace,
                   k_source=source_charge,k_obsW=obsW,k_eta_p=eta_p,
                   k_density_matrices=rho,k_e=e,k_U=U)
    return res


def response(q,theta,mu,p=Model(),eta_p0=.013164239138,
             eta_phi=.013164239138,alpha=0.,step_bar=1e-6,
             spin=IX,beta=0,obs=PA@IZ):
    r1=moments(q,theta,mu,+step_bar,p,spin,beta,obs,eta_p0,eta_phi,alpha)
    r0=moments(q,theta,mu,-step_bar,p,spin,beta,obs,eta_p0,eta_phi,alpha)
    pref=-p.g*MU_B/(2*step_bar)
    fields=('rta','bgk','correction','delta_mu','n_rta','n_bgk','source_total')
    ans={f'gamma_{k}':pref*(r1[k]-r0[k]) for k in fields}
    ans.update({'plus':r1,'minus':r0})
    return ans

if __name__=='__main__':
    for n in (101,251,401):
        q=points(n)
        for a in (0,.4):
            R=response(q,.15,.03,alpha=a)
            print(n,a,R['gamma_rta'],R['gamma_bgk'],R['gamma_correction'],R['gamma_delta_mu'],flush=True)
