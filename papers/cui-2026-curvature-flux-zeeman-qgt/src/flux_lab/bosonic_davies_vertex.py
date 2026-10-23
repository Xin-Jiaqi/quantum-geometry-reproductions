"""Round31: finite interacting Fock system coupled to declared *bosonic* spectrum.

H_SB=A \otimes B, B=sum_q g_q(b_q+b_q^dag), with J(w)=
g2*w*exp(-w/cutoff), w>0. Thermal oscillator correlations yield
Gamma(Delta)=2pi J(abs(Delta)) [n_B if Delta>0, n_B+1 if Delta<0].
A number-conserving Fock density operator creates within-N Davies/secular
transitions; zero-frequency collective dephasing uses Gamma(0)=2pi*g2/beta.

This is a microscopic-*model-derived* weak-coupling Davies dissipator under
Born, Markov, secular and continuum-spectral assumptions. Not an ab-initio
phonon spectrum, diagrammatic self-energy/vertex or full Cui metal KME.
No Lamb shift (a distinct principal-value bath renormalization) is claimed.
"""
from __future__ import annotations
import numpy as np
from .fermionic_fock_bath import (fock_soc_hubbard,grand_gibbs,
                                 unitary_phase,bath_contact)
from .thermal_lindblad_ward import dissipator,liouvillian

def bose_occupation(w,beta):
    w=float(w);beta=float(beta)
    if not np.isfinite(w) or w<=0 or not np.isfinite(beta) or beta<=0:
        raise ValueError("strictly positive boson energy and inverse temperature")
    return float(1/np.expm1(np.clip(beta*w,0,700)))

def ohmic_spectral_density(w,*,g2=.045,cutoff=2.4):
    w=float(w)
    if not np.isfinite(w) or w<0 or not np.isfinite(g2) or g2<0 or not np.isfinite(cutoff) or cutoff<=0:
        raise ValueError("nonnegative frequency/g2 and positive finite cutoff required")
    return g2*w*np.exp(-w/cutoff)

def oscillator_fourier_rate(delta,*,beta=.83,g2=.045,cutoff=2.4):
    """Delta=E_final-E_initial; + absorption, - emission."""
    delta=float(delta)
    if not np.isfinite(delta) or not np.isfinite(beta) or beta<=0:
        raise ValueError("finite Delta and positive beta")
    # Ohmic zero-frequency limiting KMS spectral weight.
    if delta==0:return 2*np.pi*g2/beta if g2>=0 and cutoff>0 and np.isfinite(g2) and np.isfinite(cutoff) else _bad_spectrum()
    w=abs(delta)
    J=ohmic_spectral_density(w,g2=g2,cutoff=cutoff)
    n=bose_occupation(w,beta)
    return 2*np.pi*J*(n if delta>0 else n+1)

def _bad_spectrum():
    raise ValueError("invalid bosonic spectral parameters")

def finite_oscillator_correlator(w,*,beta=.83,levels=36):
    """Independent actual finite Bose-oscillator matrix trace:
    C(t)=<B(t)B(0)>=emission*exp(-iwt)+absorption*exp(+iwt)
    for B=b+b^dag; the thermal Fourier ratio is KMS.
    """
    if not isinstance(levels,int) or levels<3:raise ValueError("oscillator truncation >=3")
    n=bose_occupation(w,beta)
    inds=np.arange(levels)
    weights=np.exp(-beta*w*inds)
    weights/=weights.sum()
    b=np.zeros((levels,levels),complex)
    for j in range(1,levels):b[j-1,j]=np.sqrt(j)
    B=b+b.conj().T
    hs=np.diag(inds*w)
    dm=np.diag(weights)
    # Independently sum exact Fourier coefficients via matrix elements.
    # emission: (i,j) transition energy E_j-E_i=+w, p_i |B_ij|^2
    absorption=0.;emission=0.
    for i in range(levels):
        for j in range(levels):
            if j==i+1:emission+=weights[i]*abs(B[i,j])**2
            if j==i-1:absorption+=weights[i]*abs(B[i,j])**2
    return dict(emission=emission,absorption=absorption,
                occupation=n,b=b,H=hs,rho=dm)

def _validate_fock(H,N,A,beta,g2,cutoff):
    H=np.asarray(H,complex);N=np.asarray(N,complex);A=np.asarray(A,complex)
    if H.shape!=(16,16) or N.shape!=H.shape or A.shape!=H.shape:
        raise ValueError("four fermion modes (16 Fock states) required")
    for x in (H,N,A):
        if not np.all(np.isfinite(x)) or not np.allclose(x,x.conj().T,atol=1e-11,rtol=0):
            raise ValueError("finite Hermitian Fock H,N,A")
    if np.max(abs(H@N-N@H))>1e-10 or np.max(abs(A@N-N@A))>1e-10:
        raise ValueError("number-preserving H and bath coupling A required")
    ohmic_spectral_density(1.0,g2=g2,cutoff=cutoff)
    if not np.isfinite(beta) or beta<=0:raise ValueError("positive beta required")
    d=np.diag(N)
    if (np.max(abs(N-np.diag(d)))>1e-12 or
        not np.allclose(d.real,np.rint(d.real),atol=1e-12) or
        np.max(abs(d.imag))>1e-12):
        raise ValueError("canonical Fock number occupation basis")
    return H,N,A

def bose_davies_fock_jumps(H,N,A,*,beta=.83,g2=.045,cutoff=2.4):
    """Strict within-total-N secular bath; zero-frequency collective jump."""
    H,N,A=_validate_fock(H,N,A,beta,g2,cutoff)
    nums=np.rint(np.diag(N).real).astype(int)
    es=[];ns=[];vectors=[]
    for n in range(5):
        ids=np.flatnonzero(nums==n)
        if len(ids)!=[1,4,6,4,1][n]:raise ValueError("unexpected fermion sector count")
        e,u=np.linalg.eigh(H[np.ix_(ids,ids)])
        if len(e)>1 and np.min(np.diff(e))<1e-9:
            raise ValueError("degenerate Bohr manifold requires invariant block secularization")
        for j in range(len(e)):
            v=np.zeros(16,complex);v[ids]=u[:,j]
            es.append(e[j]);ns.append(n);vectors.append(v)
    e=np.array(es);ns=np.array(ns);U=np.column_stack(vectors)
    Ab=U.conj().T@A@U
    W=np.zeros((16,16));jumps=[];channels=[]
    for m in range(16):
        for n in range(16):
            if m==n or ns[m]!=ns[n]:continue
            r=oscillator_fourier_rate(e[m]-e[n],beta=beta,g2=g2,cutoff=cutoff)*abs(Ab[m,n])**2
            W[m,n]=r
            if r>1e-16:
                jumps.append(np.sqrt(r)*np.outer(U[:,m],U[:,n].conj()))
                channels.append((m,n))
    # A(omega=0), a single Hermitian collective diagonal bath jump.
    diag=np.real(np.diag(Ab))
    Lzero=np.sqrt(oscillator_fourier_rate(0,beta=beta,g2=g2,cutoff=cutoff))*(
         (U*diag)@U.conj().T)
    if np.max(abs(Lzero))>1e-14:jumps.append(Lzero)
    if not jumps:raise ValueError("no physical bath coupling channels")
    return dict(energies=e,numbers=ns,eigenvectors=U,
                rates=W,jumps=np.asarray(jumps),transition_count=len(channels),
                zero_frequency_jump=Lzero,beta=beta,g2=g2,cutoff=cutoff,
                born_markov_secular=True,lamb_shift_included=False)

def bath_magnetic_derivative(data,*,beta=.83,g2=.045,cutoff=2.4,
                             eps=3e-5,mu=.16):
    """Centered finite Zeeman derivative of bath with independently diagonalized H+/-B."""
    if eps<=0 or not np.isfinite(eps):raise ValueError("positive finite field step required")
    H=data['H'];N=data['N'];A=data['coupling'];V=data['zeeman']
    rho=grand_gibbs(H,N,beta=beta,mu=mu)
    def shifted(sign):
        hp=H+sign*eps*V
        rb=grand_gibbs(hp,N,beta=beta,mu=mu)
        jb=bose_davies_fock_jumps(hp,N,A,beta=beta,g2=g2,cutoff=cutoff)['jumps']
        return rb,jb
    plus,Jp=shifted(1);minus,Jm=shifted(-1)
    base=bose_davies_fock_jumps(H,N,A,beta=beta,g2=g2,cutoff=cutoff)['jumps']
    dr=(plus-minus)/(2*eps)
    contact=(dissipator(rho,Jp)-dissipator(rho,Jm))/(2*eps)
    return dict(dr_equilibrium=dr,collision_contact=contact,
                equilibrium_collision_state=dissipator(dr,base),
                residual=dissipator(dr,base)+contact)

def linear_response(data,*,z=.38+.09j,mode='magnetic',beta=.83,mu=.16,
                    g2=.045,cutoff=2.4,eps=3e-5,include_bath_contact=True):
    """Retarded finite-z Liouville response, exact N² solve in Fock basis.

    mode magnetic: finite-B derived D'_B(rho) included in source, V=physical spin.
    mode gauge: pure EM G at finite z, with matching scalar+Peierls and
    bath-gauge contact; expected dr=i[G,rho], charge and physical current zero.
    """
    z=complex(z)
    if not np.isfinite(z) or z.imag<=0:raise ValueError("retarded Im z>0 required")
    H=data['H'];N=data['N'];G=data['gauge']
    rho=grand_gibbs(H,N,beta=beta,mu=mu)
    spec=bose_davies_fock_jumps(H,N,data['coupling'],beta=beta,g2=g2,cutoff=cutoff)
    J=spec['jumps']
    if mode=='gauge':
        V=1j*(G@H-H@G)+1j*z*G
        contact=bath_contact(rho,J,G)
    elif mode=='magnetic':
        V=data['zeeman']
        contact=bath_magnetic_derivative(data,beta=beta,g2=g2,cutoff=cutoff,eps=eps,mu=mu)['collision_contact']
    else:raise ValueError("mode gauge or magnetic required")
    source=-1j*(V@rho-rho@V)
    if include_bath_contact:source+=contact
    n=len(H);Id=np.eye(n*n)
    L=liouvillian(H,J)
    # Exactly (-i z - L)dr = -i[V,rho] + D'[rho].
    dr=np.linalg.solve(-1j*z*Id-L,source.reshape(n*n,order="F")).reshape((n,n),order="F")
    coherent_current=-1j*(H@data['left']-data['left']@H)
    delta_J=(1j*(G@coherent_current-coherent_current@G)
             if mode=='gauge' else np.zeros_like(H))
    current=np.trace(coherent_current@dr)+np.trace(delta_J@rho)
    density=np.array([np.trace(P@dr) for P in (data['left'],data['right'])])
    bath_source=np.array([np.trace(P@(dissipator(dr,J)+(contact if include_bath_contact else 0)))
                         for P in (data['left'],data['right'])])
    continuity=np.array([-1j*z*density[0]+current-bath_source[0],
                         -1j*z*density[1]-current-bath_source[1]])
    return dict(dr=dr,density=density,current=current,ward=continuity,
                bath_source=bath_source,source=source,contact=contact,
                rho=rho,jumps=J,z=z,mode=mode)
