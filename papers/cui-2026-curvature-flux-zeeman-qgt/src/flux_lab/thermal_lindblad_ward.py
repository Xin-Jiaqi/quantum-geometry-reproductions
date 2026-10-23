"""Round29: bounded CP/Gibbs/charge Ward contrasts for explicit finite GKSL baths.

This is an independently specified TWO-SITE spin-1/2 toy, not a transport
derivation or ab-initio bath for Cui. A canonical *one-particle* Gibbs state,
not fermionic Fermi-Dirac occupations, is used in the thermal comparison.

Bath A: local P_i jumps, CP and each site charge conserving, but not Gibbs
preserving if [P_i,H] != 0. Bath B: rank-one energy transitions with pairwise
canonical detailed balance, CP and Gibbs preserving, but bath transfers
population between sites and requires an explicit local bath source.

Physical local U(1) covariance of Bath B assumes that ALL jump operators
transform along with H. Freezing the jumps breaks the bath gauge vertex.
No frequency-dependent microscopic conserving KME or SI normalization.
"""
from __future__ import annotations
import numpy as np

I2=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)


def two_site_soc():
    """Explicit spin-orbit 4-level H, onsite site projectors, bath coupling."""
    H=np.zeros((4,4),complex)
    H[:2,:2]=.43*SZ+.11*SX
    H[2:,2:]=-.29*SZ+.13*SY+.17*I2
    T=.63*I2+.19j*SY+.12*SZ
    H[:2,2:]=T;H[2:,:2]=T.conj().T
    P0=np.diag([1,1,0,0]).astype(complex)
    P1=np.eye(4,dtype=complex)-P0
    A=P0-P1+.21*np.kron(I2,SX)+.13*np.kron(I2,SZ)
    return H,np.stack((P0,P1)),A


def canonical_gibbs(H,beta=.8):
    H=np.asarray(H,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or not np.isfinite(beta) or beta<=0:
        raise ValueError("square H and positive finite beta required")
    if not np.all(np.isfinite(H)) or not np.allclose(H,H.conj().T,atol=1e-12):
        raise ValueError("finite Hermitian H required")
    e,U=np.linalg.eigh(H)
    w=np.exp(-beta*(e-e.min()))
    return (U*(w/w.sum()))@U.conj().T


def local_lindblad(projectors,gamma=.18):
    """Local onsite dephasing P_i, exact local-charge conservation."""
    P=np.asarray(projectors,complex)
    if P.ndim!=3 or P.shape[1]!=P.shape[2] or gamma<0 or not np.isfinite(gamma):
        raise ValueError("valid local projectors and gamma>=0 required")
    n=P.shape[1]
    if (not np.allclose(sum(P),np.eye(n),atol=1e-12) or
        any(not np.allclose(p,p.conj().T,atol=1e-12) or
            not np.allclose(p@p,p,atol=1e-12) for p in P)):
        raise ValueError("orthogonal complete site projectors required")
    return np.sqrt(gamma)*P


def thermal_transition_jumps(H,A,beta=.8,gamma=.16):
    """Constructed canonical KMS pair rates W[m,n] for n -> m.

    W[m,n] = gamma |<m|A|n>|^2 exp[-beta(E_m-E_n)/2], m!=n.
    Thus W[m,n] p[n] == W[n,m] p[m], p[n]~exp(-beta E_n).
    This *postulated* bath coupling is not an electron-phonon derivation.
    """
    H=np.asarray(H,complex);A=np.asarray(A,complex)
    if (H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2 or
        not np.all(np.isfinite(H)) or not np.allclose(H,H.conj().T,atol=1e-12,rtol=0)):
        raise ValueError("finite Hermitian nondegenerate bath Hamiltonian required")
    if A.shape!=H.shape or not np.all(np.isfinite(A)) or not np.allclose(A,A.conj().T,atol=1e-12):
        raise ValueError("finite same-basis Hermitian bath coupling required")
    if not np.isfinite(beta) or beta<=0 or not np.isfinite(gamma) or gamma<0:
        raise ValueError("positive beta, finite nonnegative gamma required")
    e,U=np.linalg.eigh(H)
    if np.min(np.diff(e))<=1e-9:
        raise ValueError("degenerate bath spectrum requires invariant cluster jumps")
    B=U.conj().T@A@U;n=len(e)
    rates=np.zeros((n,n))
    jumps=[]
    for m in range(n):
        for k in range(n):
            if m==k:continue
            rate=gamma*abs(B[m,k])**2*np.exp(-.5*beta*(e[m]-e[k]))
            rates[m,k]=rate
            jumps.append(np.sqrt(rate)*np.outer(U[:,m],U[:,k].conj()))
    return dict(jumps=np.asarray(jumps),rates=rates,energies=e,eigenvectors=U)


def dissipator(X,jumps):
    """Actual GKSL jump sum, complex-operator linear, trace preserving."""
    X=np.asarray(X,complex);J=np.asarray(jumps,complex)
    if J.ndim!=3 or J.shape[1:]!=X.shape:
        raise ValueError("jumps must be (M,N,N) in same basis as X")
    out=np.zeros_like(X)
    for L in J:
        G=L.conj().T@L
        out+=L@X@L.conj().T-.5*(G@X+X@G)
    return out


def liouvillian(H,jumps):
    """Independent Kronecker N^2 generator for column-major vectorization."""
    n=len(H);I=np.eye(n,dtype=complex)
    out=-1j*(np.kron(I,H)-np.kron(H.T,I))
    for L in jumps:
        G=L.conj().T@L
        out+=np.kron(L.conj(),L)-.5*(np.kron(I,G)+np.kron(G.T,I))
    return out


def master_rhs(H,X,jumps):
    return -1j*(H@X-X@H)+dissipator(X,jumps)


def population_rhs(rates,p):
    """Independent classical detailed-balance population master equation."""
    rates=np.asarray(rates,float);p=np.asarray(p,float)
    return rates@p-p*np.sum(rates,axis=0)


def bond_current(H):
    """Physical coherent current from site 0 to site 1 in synthetic units."""
    J=np.zeros_like(H)
    J[:2,2:]=1j*H[:2,2:]
    J[2:,:2]=-1j*H[2:,:2]
    return J


def continuity_terms(H,X,projectors,jumps):
    """d n_i/dt + div J_i = Tr(P_i D(X)); bath current NOT invented."""
    dn=np.array([np.trace(P@master_rhs(H,X,jumps)) for P in projectors])
    current=np.trace(bond_current(H)@X)
    div=np.array([current,-current])
    source=np.array([np.trace(P@dissipator(X,jumps)) for P in projectors])
    return dict(density_rate=dn,divergence=div,bath_source=source,
                local_residual=dn+div-source)


def gauge_rotate(X,G,amplitude):
    G=np.asarray(G,complex)
    if not np.allclose(G,np.diag(np.diag(G)),atol=1e-12) or not np.allclose(G,G.conj().T,atol=1e-12):
        raise ValueError("diagonal Hermitian local gauge required")
    U=np.diag(np.exp(1j*amplitude*np.diag(G)))
    return U@X@U.conj().T


def bath_gauge_vertex(X,jumps,G):
    """Exact derivative of gauge-transformed dissipator at fixed X.

    D^phi[X] = U D[U^dag X U] U^dag,
    d_phi D^phi[X]|0 = i[G,D[X]] - D[i[G,X]].
    Gauge derivative of jumps is essential, unlike frozen jump rates.
    """
    DX=dissipator(X,jumps)
    dX=1j*(G@X-X@G)
    return 1j*(G@DX-DX@G)-dissipator(dX,jumps)


def evolved_rk4(H,X,jumps,dt=.01,steps=150):
    if not np.isfinite(dt) or dt<=0 or not isinstance(steps,int) or steps<1:
        raise ValueError("positive RK4 step and integer steps required")
    X=np.asarray(X,complex).copy()
    F=lambda y:master_rhs(H,y,jumps)
    for _ in range(steps):
        k1=F(X);k2=F(X+dt*k1/2)
        k3=F(X+dt*k2/2);k4=F(X+dt*k3)
        X+=dt*(k1+2*k2+2*k3+k4)/6
    return X


def evolved_dense(H,X,jumps,time=1.5):
    """Independent exact finite-matrix exponential by eigensystem."""
    n=len(H)
    L=liouvillian(H,jumps)
    w,S=np.linalg.eig(L)
    vec=np.asarray(X,complex).reshape(n*n,order="F")
    result=(S*np.exp(w*time))@np.linalg.solve(S,vec)
    return result.reshape((n,n),order="F")
