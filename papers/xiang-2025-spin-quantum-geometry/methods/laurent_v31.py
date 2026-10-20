"""Adiabatic eta -> 0 coefficient audit for general isolated spectral clusters.

Independent NumPy migration of offline V17--V19 spectral projector derivation.
NOT a scattering-driven metallic DC intrinsic solver.

X = d_k f(H); Pi Y = sum(Pa Y Pa); Qj Y = sum(a!=b Pa Y Pb/[i(Ea-Eb)]**j).
For rho2 = L2^-1 d_k L1^-1 X with L1=eta+i[H,.], L2=2eta+i[H,.],

rho2 = A/eta**2 + B/eta + C + O(eta)
A = .5 Pi d_k(Pi X)
B = .5 Pi d_k(Q1 X) + Q1 d_k(Pi X)
C = -.5 Pi d_k(Q2 X) + Q1 d_k(Q1 X) - 2 Q2 d_k(Pi X)

The C value is scheme-specific; never identify it as a validated
physical intrinsic tau^0 coefficient without a separate DC protocol.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .finite_eta import equilibrium, fd4, hermitian, resolvent

@dataclass(frozen=True)
class SpectralCluster:
    energies: tuple[float,...]
    projectors: tuple[np.ndarray,...]
    gap_min: float

def spectral_clusters(H,atol=2e-8,min_gap=1e-5):
    """Return full isolated-energy-subspace projectors; never split doublets."""
    H=hermitian(H)
    if atol<=0 or min_gap<=atol:raise ValueError('0<atol<min_gap required')
    e,U=np.linalg.eigh(H);groups=[];start=0
    for j in range(1,len(e)):
        if e[j]-e[j-1]>atol:
            groups.append(slice(start,j));start=j
    groups.append(slice(start,len(e)))
    if len(groups)<2:raise ValueError('at least two isolated clusters required')
    gaps=[float(e[s2.start]-e[s1.stop-1]) for s1,s2 in zip(groups[:-1],groups[1:])]
    if min(gaps)<=min_gap:raise ValueError('intercluster gap too small')
    if any(np.ptp(e[s])>atol for s in groups):raise ValueError('cluster spread')
    Ps=tuple(U[:,s]@U[:,s].conj().T for s in groups)
    Es=tuple(float(e[s].mean()) for s in groups)
    return SpectralCluster(Es,Ps,float(min(gaps)))

def block(Y,clusters,power=0,diagonal=False):
    if power not in (0,1,2):raise ValueError('power must be 0,1,2')
    out=np.zeros_like(Y,dtype=complex)
    for a,(ea,Pa) in enumerate(zip(clusters.energies,clusters.projectors)):
        for b,(eb,Pb) in enumerate(zip(clusters.energies,clusters.projectors)):
            if (a==b)!=diagonal:continue
            if power and a==b:continue
            out+=Pa@Y@Pb/(1j*(ea-eb))**power if power else Pa@Y@Pb
    return out

def _at(Hfun,kx,ky,mu,T,axis,step,atol,min_gap):
    k=np.array([kx,ky],float)
    def at(z):
        q=k.copy();q[axis]=z
        return hermitian(Hfun(*q))
    c=spectral_clusters(at(k[axis]),atol,min_gap)
    X=fd4(lambda z:equilibrium(at(z),mu,T),k[axis],step)
    return c,X,block(X,c,diagonal=True),block(X,c,power=1),block(X,c,power=2)

def adiabatic_laurent(Hfun,kx,ky,mu=0.,temperature=.2,axis=0,
                      step=4e-4,atol=2e-8,min_gap=1e-5,outer_rate=2.):
    """Extract A,B,C in orbitals with fixed ranks along the k derivative stencil."""
    if axis not in (0,1):raise ValueError('axis 0 or 1')
    if outer_rate<=0:raise ValueError('outer_rate positive')
    if step<=0 or temperature<=0:raise ValueError('positive step and temperature required')
    k=np.array([kx,ky],float)
    def parts(z):
        q=k.copy();q[axis]=z
        return _at(Hfun,*q,mu,temperature,axis,step,atol,min_gap)
    c,X,pop,Y,Z=parts(k[axis])
    ranks=tuple(round(np.trace(P).real) for P in c.projectors)
    def partial(z,idx):
        d=parts(z)
        if tuple(round(np.trace(P).real) for P in d[0].projectors)!=ranks:
            raise ValueError('cluster ranks changed under k-derivative stencil')
        return d[idx]
    dp=fd4(lambda z:partial(z,2),k[axis],step)
    dy=fd4(lambda z:partial(z,3),k[axis],step)
    dz=fd4(lambda z:partial(z,4),k[axis],step)
    diagonal=-block(dz,c,diagonal=True)/outer_rate
    mixing=block(dy,c,power=1)
    occupation=-outer_rate*block(dp,c,power=2)
    A=block(dp,c,diagonal=True)/outer_rate
    B=block(dy,c,diagonal=True)/outer_rate+block(dp,c,power=1)
    C=diagonal+mixing+occupation
    return {'A':A,'B':B,'C_adiabatic':C,'diagonal':diagonal,
            'mixing':mixing,'occupation':occupation,
            'delta_printed_candidate':occupation/2 if outer_rate==2. else None,
            'C_half_population_candidate':diagonal+mixing+occupation/2 if outer_rate==2. else None,
            'cluster_ranks':ranks,'gap_min':c.gap_min,'outer_rate':float(outer_rate),
            'protocol':f'adiabatic_L1_eta_L2_{outer_rate:g}eta_zero_frequency',
            'scientific_scope':'finite eta Laurent constant is NOT established physical intrinsic DC'}

def independent_population_delta(Hfun,kx,ky,mu=0.,temperature=.2,
                                 axis=0,step=4e-4,atol=2e-8,min_gap=1e-5):
    """Compute -Q2 d_k(Pi X) via dH/dk and f'(Ea), not derivatives of Pi X."""
    if axis not in (0,1):raise ValueError('invalid axis')
    k=np.array([kx,ky],float)
    def at(z):
        q=k.copy();q[axis]=z
        return hermitian(Hfun(*q))
    H=at(k[axis]);c=spectral_clusters(H,atol,min_gap)
    Hp=fd4(at,k[axis],step)
    fp=[]
    for ea,Pa in zip(c.energies,c.projectors):
        ea_prime=float(np.trace(Pa@Hp).real/round(np.trace(Pa).real))
        f=1/(1+np.exp(np.clip((ea-mu)/temperature,-700,700)))
        fp.append(-f*(1-f)*ea_prime/temperature)
    out=np.zeros_like(H)
    for a,(ea,Pa) in enumerate(zip(c.energies,c.projectors)):
        for b,(eb,Pb) in enumerate(zip(c.energies,c.projectors)):
            if a==b:continue
            dPi=(fp[b]-fp[a])*(Pa@Hp@Pb)/(eb-ea)
            out-=dPi/(1j*(ea-eb))**2
    return out

def reconstruct_response(coefficients,eta,observable=None):
    """Truncated Laurent, not an intrinsic tau^0 extraction."""
    if eta<=0:raise ValueError('eta positive')
    matrix=coefficients['A']/eta**2+coefficients['B']/eta+coefficients['C_adiabatic']
    if observable is None:return matrix
    O=hermitian(observable)
    return float(np.trace(O@matrix).real)

def explicit_three_cluster_term(Hfun,kx,ky,mu=0.,temperature=.2,axis=0,
                                step=4e-4,atol=2e-8,min_gap=1e-5):
    """Intermediate c!=a,b clusters in Q1 d_k(Q1X), with projectors only."""
    k=np.array([kx,ky],float)
    if axis not in (0,1):raise ValueError('invalid axis')
    def at(z):
        q=k.copy();q[axis]=z
        return hermitian(Hfun(*q))
    H=at(k[axis]);c=spectral_clusters(H,atol,min_gap)
    X=fd4(lambda z:equilibrium(at(z),mu,temperature),k[axis],step)
    Y=block(X,c,power=1);Hp=fd4(at,k[axis],step)
    result=np.zeros_like(H)
    for a,(ea,Pa) in enumerate(zip(c.energies,c.projectors)):
        for b,(eb,Pb) in enumerate(zip(c.energies,c.projectors)):
            if a==b:continue
            between=np.zeros_like(H)
            for i,(ec,Pc) in enumerate(zip(c.energies,c.projectors)):
                if i in (a,b):continue
                rac=-1j*(Pa@Hp@Pc)/(ea-ec)
                rcb=-1j*(Pc@Hp@Pb)/(ec-eb)
                between+=rac@Pc@Y@Pb-Pa@Y@Pc@rcb
            result+=-1j*between/(1j*(ea-eb))
    return result

def second_order_outer_rate(Hfun,kx,ky,eta,outer_rate=2.,mu=0.,
                            temperature=.2,axis=0,step=4e-4):
    """Exact finite-eta check; outer_rate!=2 is COUNTERFACTUAL only."""
    if eta<=0 or outer_rate<=0 or axis not in (0,1):raise ValueError('bad eta/rate/axis')
    k=np.array([kx,ky],float)
    def at(x):
        q=k.copy();q[axis]=x
        return hermitian(Hfun(*q))
    def first(x):
        H=at(x)
        X=fd4(lambda y:equilibrium(at(y),mu,temperature),x,step)
        return -resolvent(H,X,eta=eta,order=1,omega=0.)
    drho1=fd4(first,k[axis],step)
    E,U=np.linalg.eigh(at(k[axis]))
    out=U@((U.conj().T@drho1@U)/(outer_rate*eta+1j*(E[:,None]-E[None,:])))@U.conj().T
    return -out
