"""Local NONDEGENERATE interband-only GHE/KME protocol for Cui et al. (2026).

H, D_kH and three *physical* dimensionless Pauli-spin operators must be in
the same fixed orthonormal orbital basis. E drive +E.r, magnetic drive
H_B=-kappa B.sigma; J=D_kH, M=kappa sigma, hbar=1 (natural units).
Source Fourier omega_source = -omega_drive; never silently switch signs.
No intraband, non-Abelian, collision, bulk-position, or SI claim.
"""
from __future__ import annotations
import numpy as np


def _fermi(e, mu, temperature):
    if temperature < 0: raise ValueError("temperature must be nonnegative")
    if temperature == 0: return (np.asarray(e) < mu).astype(float)
    return 1/(1+np.exp(np.clip((np.asarray(e)-mu)/temperature,-700,700)))


def _data(H, dH, spin, *, gap_tolerance=1e-9):
    H=np.asarray(H,complex); dH=np.asarray(dH,complex); spin=np.asarray(spin,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1] or H.shape[0]<2:
        raise ValueError("H must be a square N>=2 matrix")
    N=H.shape[0]
    if dH.shape!=(2,N,N) or spin.shape!=(3,N,N):
        raise ValueError("Require two D_kH and three physical-spin matrices")
    for A in [H,*dH,*spin]:
        if not np.all(np.isfinite(A)): raise ValueError("Nonfinite input")
        if not np.allclose(A,A.conj().T,rtol=0,atol=1e-10):
            raise ValueError("All operators must be Hermitian in the same basis")
    e,U=np.linalg.eigh(H)
    if np.min(np.diff(e))<=gap_tolerance:
        raise ValueError("Degenerate spectrum: non-Abelian protocol required")
    v=np.einsum("ni,anm,mj->aij",U.conj(),dH,U,optimize=True)
    s=np.einsum("ni,bnm,mj->bij",U.conj(),spin,U,optimize=True)
    r=np.zeros_like(v)
    for n in range(N):
        for m in range(N):
            if n!=m:r[:,n,m]=v[:,n,m]/(1j*(e[n]-e[m]))
    if not np.allclose(r,np.swapaxes(r.conj(),-1,-2),rtol=0,atol=1e-10):
        raise AssertionError("Interband position must be Hermitian")
    R=np.einsum("in,anm,jm->aij",U,r,U.conj(),optimize=True)
    return dict(e=e,U=U,v=v,s=s,r=r,R=R)


def interband_kubo(H,dH,spin,*,mu=0.,temperature=.025,z=.04,kappa=1.):
    """Directed band pairs; outputs chi_JB(2,3), alpha_ME(3,2)."""
    D=_data(H,dH,spin);e,v,s,r=(D[k] for k in ("e","v","s","r"))
    f=_fermi(e,mu,temperature);z=complex(z)
    JB=np.zeros((2,3),complex);ME=np.zeros((3,2),complex)
    for n in range(len(e)):
        for m in range(len(e)):
            if n==m:continue
            q=f[n]-f[m];den=z-(e[m]-e[n])
            if abs(den)<1e-10:
                if abs(q)>1e-10:raise ValueError("On-shell pole")
                continue
            JB-=kappa*q*np.outer(v[:,n,m],s[:,m,n])/den
            ME+=kappa*q*np.outer(s[:,n,m],r[:,m,n])/den
    return JB,ME


def interband_liouville(H,dH,spin,*,mu=0.,temperature=.025,z=.04,kappa=1.):
    """Independent orbital-basis N² Liouvillian (no spectral response sum)."""
    D=_data(H,dH,spin)
    H=np.asarray(H,complex);dH=np.asarray(dH,complex);spin=np.asarray(spin,complex)
    e,U,R=D["e"],D["U"],D["R"];z=complex(z);N=len(e)
    if abs(z)<1e-12:raise ValueError("Full Liouvillian singular at z=0")
    f=_fermi(e,mu,temperature);rho0=(U*f)@U.conj().T
    L=np.kron(np.eye(N),H)-np.kron(H.T,np.eye(N))
    M=z*np.eye(N*N,dtype=complex)-L
    JB=np.zeros((2,3),complex);ME=np.zeros((3,2),complex)
    for b,S in enumerate(spin):
        Q=-kappa*(S@rho0-rho0@S)
        d=np.linalg.solve(M,Q.reshape(-1,order="F")).reshape((N,N),order="F")
        for a,V in enumerate(dH):JB[a,b]=np.trace(V@d)
    for a,P in enumerate(R):
        Q=P@rho0-rho0@P
        d=np.linalg.solve(M,Q.reshape(-1,order="F")).reshape((N,N),order="F")
        for b,S in enumerate(spin):ME[b,a]=kappa*np.trace(S@d)
    return JB,ME


def static_subtracted_gyrotropic(H,dH,spin,*,mu=0.,temperature=.025,
                                 omega_drive=.04,eta=0.,kappa=1.):
    """Subtract chi_JB(i eta) at the SAME causal regulator."""
    if eta<0:raise ValueError("eta must be nonnegative")
    opts=dict(mu=mu,temperature=temperature,kappa=kappa)
    return (interband_kubo(H,dH,spin,z=complex(omega_drive,eta),**opts)[0]
            -interband_kubo(H,dH,spin,z=1j*eta,**opts)[0])


def reciprocal_residual(H,dH,spin,*,mu=0.,temperature=.025,
                        omega_drive=.04,eta=0.,kappa=1.):
    """Exact local interband identity, including static alpha_ME.

    Delta chi_JB(z)=i[z alpha_ME(-z)-z0 alpha_ME(-z0)]^T,
    z=omega_drive+i eta, z0=i eta. For eta>0, alpha(-z) is advanced;
    replacing it with retarded alpha(-omega+i eta) is invalid.
    """
    z=complex(omega_drive,eta);z0=1j*eta
    opts=dict(mu=mu,temperature=temperature,kappa=kappa)
    lhs=static_subtracted_gyrotropic(H,dH,spin,omega_drive=omega_drive,eta=eta,**opts)
    a=interband_kubo(H,dH,spin,z=-z,**opts)[1]
    a0=interband_kubo(H,dH,spin,z=-z0,**opts)[1]
    rhs=1j*(z*a-z0*a0).T
    wrong=1j*(z*(a-a0)).T
    return dict(max_residual=float(np.max(abs(lhs-rhs))),
                max_wrong_double_subtraction=float(np.max(abs(lhs-wrong))),
                delta_chi_JB=lhs,rhs=rhs,alpha_full=a,alpha_static=a0,z=z,z0=z0)


def source_four_sector_expansion(H,dH,spin,*,mu=0.,temperature=.025,
                                 omega_source=.01,kappa=1.):
    """Eq.(16) source Fourier frequency through O(omega_source²).
    Derived directly from T_nm = r_nm sigma_mn and eps_nm=e_n-e_m.
    """
    D=_data(H,dH,spin);e,r,s=(D[k] for k in ("e","r","s"))
    f=_fermi(e,mu,temperature)
    c1=np.zeros((2,3),complex);c2=np.zeros((2,3),complex)
    for n in range(len(e)):
        for m in range(len(e)):
            if n==m:continue
            eps=e[n]-e[m]
            T=np.outer(r[:,n,m],s[:,m,n])
            c1+=-2j*kappa*f[n]*T.real/eps
            c2+=2*kappa*f[n]*T.imag/eps**2
    return dict(sigma_source_omega=omega_source*c1+omega_source**2*c2,
                linear_coefficient=c1,quadratic_coefficient=c2,
                source_frequency=omega_source,drive_frequency=-omega_source)


def source_four_sector_geometric(H,dH,spin,*,mu=0.,temperature=.025,
                                 omega_source=.01,kappa=1.):
    """Independent explicit Eq.(16) gN,gA,OmegaN,OmegaA recombination.

    Embed a 2D position r in xyz with r_z=0, yielding a 2x3 output.
    The individual sector names describe algebraic pieces, not
    experimentally separable observables without symmetry analysis.
    """
    D=_data(H,dH,spin);e,r,s=(D[k] for k in ("e","r","s"))
    f=_fermi(e,mu,temperature)
    out={k:np.zeros((2,3),complex) for k in ("gN","gA","omegaN","omegaA")}
    for n in range(len(e)):
        for m in range(len(e)):
            if n==m:continue
            eps=e[n]-e[m]
            T=np.zeros((3,3),complex)
            T[:2,:]=np.outer(r[:,n,m],s[:,m,n])
            gn=((T+T.T)/2).real
            ga=((T+T.T)/2).imag
            crossA=-(T-T.T).real
            crossN=-(T-T.T).imag
            out["omegaA"]+=kappa*f[n]*(1j*omega_source/eps)*crossA[:2,:]
            out["gN"]+=kappa*f[n]*(-2j*omega_source/eps)*gn[:2,:]
            out["omegaN"]+=kappa*f[n]*(-omega_source**2/eps**2)*crossN[:2,:]
            out["gA"]+=kappa*f[n]*(2*omega_source**2/eps**2)*ga[:2,:]
    return {**out,"sigma_source_omega":sum(out.values()),
            "source_frequency":omega_source,"drive_frequency":-omega_source}
