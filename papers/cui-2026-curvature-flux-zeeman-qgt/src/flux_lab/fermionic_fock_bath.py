"""Round30: finite interacting four-mode FERMIONIC Fock-space thermal bath.

CAR-exact 2-site x physical spin-1/2, 16 Fock configurations and Hubbard
interaction. A chosen number-conserving Hermitian density coupling is
energy-sector decomposed into canonical pairwise KMS jumps. CP GKSL preserves
TOTAL particle number and interacting GRAND canonical rho~exp[-beta(H-mu N)].
Individual *site* charge changes under the nonlocal bath and its local
continuity law needs an explicit source. Bath gauge and magnetic contacts
must be transformed/differentiated consistently with thermal jumps.

No microscopic phonon spectrum, Pauli-blocked Boltzmann collision integral,
thermodynamic metallic limit or original Cui Eq.(19) / SI KME.
"""
from __future__ import annotations
import numpy as np
from .thermal_lindblad_ward import two_site_soc, dissipator, liouvillian

def fermion_ops(nmode=4):
    if nmode!=4:raise ValueError("bounded four-mode CAR protocol")
    dim=1<<nmode
    annihilate=np.zeros((nmode,dim,dim),complex)
    for mode in range(nmode):
        for state in range(dim):
            if (state>>mode)&1:
                target=state^(1<<mode)
                sign=(-1)**((state&((1<<mode)-1)).bit_count())
                annihilate[mode,target,state]=sign
    create=np.swapaxes(annihilate.conj(),-1,-2)
    return annihilate,create

def lift_one_body(h,annihilate=None,create=None):
    h=np.asarray(h,complex)
    if h.shape!=(4,4) or not np.all(np.isfinite(h)) or not np.allclose(h,h.conj().T,atol=1e-12):
        raise ValueError("Hermitian physical spinful four-mode one-body matrix")
    if annihilate is None:annihilate,create=fermion_ops()
    result=np.zeros((16,16),complex)
    for i in range(4):
        for j in range(4):result+=h[i,j]*create[i]@annihilate[j]
    return result

def fock_soc_hubbard(*,interaction=.48,cross=.19):
    if not np.isfinite(interaction) or not np.isfinite(cross):
        raise ValueError("finite interaction required")
    h,_,_=two_site_soc()
    c,cd=fermion_ops()
    occ=np.array([cd[i]@c[i] for i in range(4)])
    N=occ.sum(axis=0)
    left=occ[0]+occ[1];right=occ[2]+occ[3]
    H=lift_one_body(h,c,cd)
    H+=interaction*(occ[0]@occ[1]+.87*occ[2]@occ[3])
    H+=cross*left@right
    A=left-right+.29*(occ[0]-occ[1])
    spin_y=np.array([[0,-1j],[1j,0]],complex)
    H_B=lift_one_body(-.7*np.kron(np.eye(2),spin_y),c,cd)
    G=.43*left-.32*right
    return dict(H=H,N=N,left=left,right=right,occupations=occ,
                coupling=A,zeeman=H_B,gauge=G,annihilate=c,create=cd,
                interaction=interaction,cross=cross,h=h)

def grand_gibbs(H,N,*,beta=.83,mu=.16):
    H=np.asarray(H,complex);N=np.asarray(N,complex)
    if H.shape!=(16,16) or N.shape!=H.shape or not np.isfinite(beta) or beta<=0 or not np.isfinite(mu):
        raise ValueError("finite 16-state grand-canonical beta>0, mu")
    if not np.allclose(H,H.conj().T,atol=1e-11) or not np.allclose(N,N.conj().T,atol=1e-11):
        raise ValueError("Hermitian H,N")
    if np.max(abs(H@N-N@H))>1e-10:
        raise ValueError("total N must be conserved")
    e,U=np.linalg.eigh(H-mu*N)
    p=np.exp(-beta*(e-e.min()));p/=p.sum()
    return (U*p)@U.conj().T

def thermal_sector_jumps(H,N,A,*,beta=.83,gamma=.14):
    """Construct exact Fock-sector rank-one jumps from a declared bath coupling.

    W[m,n]=gamma |A_mn|^2 exp(-beta*(E_m-E_n)/2), only if Nm=Nn.
    These rates encode a *postulated* bath KMS ratio, not a microscopic
    derived electron-phonon spectral density or self-energy.
    """
    H=np.asarray(H,complex);N=np.asarray(N,complex);A=np.asarray(A,complex)
    if H.shape!=(16,16) or N.shape!=H.shape or A.shape!=H.shape:
        raise ValueError("four-mode 16-state operators required")
    for X in (H,N,A):
        if not np.all(np.isfinite(X)) or not np.allclose(X,X.conj().T,atol=1e-11):
            raise ValueError("finite Hermitian operators required")
    if (np.max(abs(H@N-N@H))>1e-10 or np.max(abs(A@N-N@A))>1e-10):
        raise ValueError("H and bath coupling must conserve total N")
    if not np.isfinite(beta) or beta<=0 or not np.isfinite(gamma) or gamma<0:
        raise ValueError("beta>0 and gamma>=0")
    dN=np.diag(N).real
    if not np.allclose(N,np.diag(np.diag(N)),atol=1e-12) or not np.allclose(dN,np.rint(dN),atol=1e-12):
        raise ValueError("canonical Fock occupation-number basis required")
    energies=[];numbers=[];cols=[]
    for sector in range(5):
        ids=np.where(np.rint(dN).astype(int)==sector)[0]
        if len(ids)!=[1,4,6,4,1][sector]:raise ValueError("wrong Fock sector dimensions")
        es,Us=np.linalg.eigh(H[np.ix_(ids,ids)])
        if len(es)>1 and np.min(np.diff(es))<1e-9:
            raise ValueError("exactly degenerate sector needs invariant cluster jumps")
        for j,e in enumerate(es):
            v=np.zeros(16,complex);v[ids]=Us[:,j]
            energies.append(e);numbers.append(sector);cols.append(v)
    E=np.asarray(energies);Ns=np.asarray(numbers)
    U=np.column_stack(cols);Ba=U.conj().T@A@U
    rates=np.zeros((16,16))
    jumps=[]
    for m in range(16):
        for n in range(16):
            if m==n or Ns[m]!=Ns[n]:continue
            r=gamma*abs(Ba[m,n])**2*np.exp(-beta*(E[m]-E[n])/2)
            rates[m,n]=r
            if r>1e-16:jumps.append(np.sqrt(r)*np.outer(U[:,m],U[:,n].conj()))
    if not jumps:raise ValueError("no thermalizing number-conserving channels")
    return dict(jumps=np.array(jumps),rates=rates,energies=E,
                numbers=Ns,eigenvectors=U)

def local_charge_balance(data,X,jumps):
    """Exact manybody site Ward: d n_left + J_coherent = bath_source_left."""
    H=data['H'];left=data['left'];right=data['right']
    X=np.asarray(X,complex)
    D=dissipator(X,jumps)
    vel=-1j*(H@X-X@H)+D
    J=-1j*(H@left-left@H)
    current=np.trace(J@X)
    dl=np.trace(left@vel);dr=np.trace(right@vel)
    src_l=np.trace(left@D);src_r=np.trace(right@D)
    return dict(current=current,density_rate=np.array([dl,dr]),
         bath_source=np.array([src_l,src_r]),
         continuity=np.array([dl+current-src_l,dr-current-src_r]),
         total_source=np.trace(data['N']@D))

def unitary_phase(X,G,amplitude):
    G=np.asarray(G,complex)
    if G.shape!=(16,16) or not np.allclose(G,np.diag(np.diag(G)),atol=1e-11) or not np.allclose(G,G.conj().T,atol=1e-11):
        raise ValueError("diagonal real Fock charge gauge")
    U=np.diag(np.exp(1j*amplitude*np.diag(G)))
    return U@np.asarray(X,complex)@U.conj().T

def bath_contact(X,jumps,G):
    """Generator-contact d_phi D_phi(X) at phi=0 with all jumps rotated."""
    D=dissipator(X,jumps)
    dX=1j*(G@X-X@G)
    return 1j*(G@D-D@G)-dissipator(dX,jumps)

def sector_gibbs_weights(E,N,*,beta=.83,mu=.16):
    a=np.exp(-beta*(np.asarray(E)-mu*np.asarray(N)-
           np.min(np.asarray(E)-mu*np.asarray(N))))
    return a/a.sum()
