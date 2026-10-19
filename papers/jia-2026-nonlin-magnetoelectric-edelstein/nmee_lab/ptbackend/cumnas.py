"""CuMnAs minimal Hamiltonian Eq (5) of Jia et al PRL 136, 226302 (2026).

Reproducibility conventions:
- q=(kx,ky) dimensionless, BZ [-pi,pi)^2 and a physical lattice constant a_nm.
- t in eV; t'=0.08t, lambda=0.8t, J=0.6t.
- H uses tau (A/B) tensor sigma (real spin Pauli).
- Energies E_- and E_+ each have EXACT double degeneracy: do not use naive nondegenerate band sums.
- Intrinsic result computed by finite-field differentiation of clean interband static
  spin-electric susceptibility with frozen occupations. This tests a full-static
  definition, NOT automatically the author's printed Eq.(1) for degenerate bands.
- Extrinsic uses the printed Eq.(2) with doublet projectors (no 0/0), with
  tau_fs and hbar in eV fs. The physical length a_nm is NOT supplied in the paper.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import expit

SX=np.array([[0,1],[1,0]],complex); SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1,-1]).astype(complex); I=np.eye(2, dtype=complex)
TX=np.kron(SX,I); TZ=np.kron(SZ,I); IX=np.kron(I,SX); IY=np.kron(I,SY); IZ=np.kron(I,SZ)
S={"x":IX,"y":IY,"z":IZ}
PA=(np.eye(4)+TZ)/2; PB=np.eye(4)-PA
MU_B=5.7883818060e-5 # eV/T
KB=8.617333262145e-5 # eV/K
HBAR_EV_FS=0.6582119569

@dataclass(frozen=True)
class Model:
    t:float=1.0
    tp_ratio:float=.08
    lam_ratio:float=.8
    j_ratio:float=.6
    g:float=10.
    a_nm:float=1.
    tau_fs:float=50.
    temperature_K:float=100.

    def params(self):return self.t,self.tp_ratio*self.t,self.lam_ratio*self.t,self.j_ratio*self.t


def points(n, offset=(0.39,0.17)):
    v=np.meshgrid(-np.pi+2*np.pi*(np.arange(n)+offset[0])/n,
                  -np.pi+2*np.pi*(np.arange(n)+offset[1])/n,indexing='ij')
    return np.stack(v,-1).reshape(-1,2)


def terms(q,theta,p=Model()):
    kx,ky=q[:,0],q[:,1];t,tp,lam,J=p.params()
    A=-2*t*np.cos(kx/2)*np.cos(ky/2)
    c=-tp*(np.cos(kx)+np.cos(ky))
    bx=J*np.cos(theta)-lam*np.sin(ky)
    by=J*np.sin(theta)+lam*np.sin(kx)
    return c,A,bx,by


def hamiltonian(q,theta,p=Model()):
    c,A,bx,by=terms(q,theta,p)
    return c[:,None,None]*np.eye(4)+A[:,None,None]*TX + bx[:,None,None]*(TZ@IX)+by[:,None,None]*(TZ@IY)


def velocity(q,theta,beta,p=Model()):
    kx,ky=q[:,0],q[:,1];t,tp,lam,J=p.params()
    if beta==0:
        dc=tp*np.sin(kx);dA=t*np.sin(kx/2)*np.cos(ky/2)
        dbx=np.zeros_like(kx); dby=lam*np.cos(kx)
    elif beta==1:
        dc=tp*np.sin(ky);dA=t*np.cos(kx/2)*np.sin(ky/2)
        dbx=-lam*np.cos(ky);dby=np.zeros_like(kx)
    else: raise ValueError('beta must be 0 or 1')
    return p.a_nm*(dc[:,None,None]*np.eye(4)+dA[:,None,None]*TX +
                dbx[:,None,None]*(TZ@IX)+dby[:,None,None]*(TZ@IY))


def energies(q,theta,p=Model()):
    c,A,bx,by=terms(q,theta,p)
    d=np.sqrt(A**2+bx**2+by**2)
    return c-d,c+d,d


def fermi(e,mu,T):
    return (e<mu).astype(float) if T==0 else expit((mu-e)/(KB*T))


def fermi_prime(e,mu,T):
    if T<=0:raise ValueError('Extrinsic fprime requires finite temperature')
    f=fermi(e,mu,T)
    return -f*(1-f)/(KB*T)


def _static_intergroup_kernel(h,v,observable,eta_eV=0.):
    """K per k: chi_interband=(f_- -f_+) K, with doublet grouping.

    K=-2 Im Σ_{n in lower,m in upper} S_nm v_mn/(e_n-e_m)^2
    Intradoublet transitions drop exactly because the unperturbed f_n coincide.
    """
    e,U=np.linalg.eigh(h)
    ud=U.conj().transpose(0,2,1)
    ve=ud@v@U; se=ud@observable@U
    den=(e[:,:2,None]-e[:,None,2:])**2+eta_eV**2
    return -2*np.imag(se[:,:2,2:]*ve[:,2:,:2].transpose(0,2,1)).__truediv__(den).sum(axis=(1,2))


def intrinsic_static_kernels(q,theta,p=Model(),Bstep_T=.02,eta_eV=0.):
    """Return k kernels for A, B and total, for local sigma_z, B_x and E_x.

    All are generated from same Bloch Hamiltonian and fixed zero-field f_+/-.
    Multiply by (f_- - f_+) and integrate 1/(a_nm^2 * Nk).
    """
    h=hamiltonian(q,theta,p); v=velocity(q,theta,0,p)
    low,up,d=energies(q,theta,p)
    if np.any(d < p.g*MU_B*Bstep_T*3):
        raise ValueError('Field step too large at nearly gapless q point; reduce Bstep or shift mesh')
    vals=[]
    for sign in (1,-1):
        hf=h-sign*Bstep_T*p.g*MU_B*IX
        vals.append(np.array([_static_intergroup_kernel(hf,v,PA@IZ,eta_eV),
                              _static_intergroup_kernel(hf,v,PB@IZ,eta_eV),
                              _static_intergroup_kernel(hf,v,IZ,eta_eV)]))
    out=(vals[0]-vals[1])/(2*Bstep_T)
    return low,up,out # (3,Nk)


def intrinsic_mu_scan(q,theta,mus,p=Model(),Bstep_T=.02,eta_eV=0.):
    low,up,k=intrinsic_static_kernels(q,theta,p,Bstep_T,eta_eV)
    occ=fermi(low[None,:],np.asarray(mus)[:,None],p.temperature_K)-fermi(up[None,:],np.asarray(mus)[:,None],p.temperature_K)
    return (occ@k.T)/(len(q)*p.a_nm**2) # n_mu,3


def _trace4(p,sa,r,sc):
    # real trace(P sa R sc), gauge invariant inside each doublet.
    return np.real(np.einsum('nij,jk,nkl,lm, nmi->n',p,sa,r,sc,np.eye(4)[None],optimize=True))


def extrinsic_fprime_weights(q,theta,p=Model(), alpha='x', gamma='y',beta=1):
    """For Eq2: 2 ∑_{m in other pair} Re(s_alpha_nm s_gamma_mn)*v_n/(e_n-e_m)
    Returns per upper/lower E_p coefficient WITHOUT occupation f'_p or global pref.
    Uses exact spectral projectors, valid despite PT-enforced spin degeneracy.
    """
    c,A,bx,by=terms(q,theta,p)
    d=np.sqrt(A*A+bx*bx+by*by)
    if np.any(d < 1.e-11): raise ValueError('Exact Dirac point: interdoublet formula undefined')
    D=A[:,None,None]*TX+bx[:,None,None]*(TZ@IX)+by[:,None,None]*(TZ@IY)
    Pplus=(np.eye(4)+D/d[:,None,None])*.5
    Pminus=(np.eye(4)-D/d[:,None,None])*.5
    vv=velocity(q,theta,beta,p)
    # Doublet-averaged group velocity: tr(P_s v)/2
    vm=np.einsum('nij,nji->n',Pminus,vv).real/2
    vp=np.einsum('nij,nji->n',Pplus,vv).real/2
    sa,sc=S[alpha],S[gamma]
    Sm=np.einsum('nij,ji->n',Pminus@sa@Pplus,sc,optimize=False).real
    Sp=np.einsum('nij,ji->n',Pplus@sa@Pminus,sc,optimize=False).real
    # Pauli spin metrics symmetric alpha/gamma for each group under trace
    wm=2*Sm*vm/(-2*d)
    wp=2*Sp*vp/(2*d)
    return c-d,c+d,wm,wp, d, Sm, Sp


def extrinsic_mu_scan(q,theta,mus,p=Model(),alpha='x',gamma='y',beta=1):
    em,ep,wm,wp,_,_,_=extrinsic_fprime_weights(q,theta,p,alpha,gamma,beta)
    mu=np.asarray(mus)[:,None]
    kernel=fermi_prime(em[None,:],mu,p.temperature_K)*wm[None,:]+fermi_prime(ep[None,:],mu,p.temperature_K)*wp[None,:]
    # g muB * tau / hbar in chosen eV, nm, fs convention
    pref=p.g*MU_B*p.tau_fs/HBAR_EV_FS
    return pref*kernel.mean(axis=1)/(p.a_nm**2)


def gapless_momenta(theta,p=Model()):
    """An exact zero of d requires A=bx=by=0. Useful only for special theta.
    """
    return None
