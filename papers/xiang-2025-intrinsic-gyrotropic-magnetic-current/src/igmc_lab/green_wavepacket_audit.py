"""Independent finite-omega checks of Xiang et al. SM III-IV.

No import of earlier Kubo, geometric pair kernels, or Fermi-contour routines.
The GF method evaluates the original SM Eq.(44) as matrix resolvents at the
spectral poles. The wavepacket method evaluates SM Eqs.(21),(24) directly at
positive and negative drive harmonics. Hamiltonian derivatives are supplied
as plain orbital-basis matrices. Frequencies are photon energies hbar*omega.
"""
from __future__ import annotations
import numpy as np

EYE=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
SPINS={'x':SX,'y':SY,'z':SZ}

def green_eq44_resolvent(H,velocity,spin,omega,mu=0.,eta=0.):
    """Return K=-J/b in model units for magnetic perturbation -b*sigma.

    Independently evaluate G1< = i sum_n f_n [Gr(e_n+omega) Sigma P_n
                               + P_n Sigma Ga(e_n-omega)] for Sigma=-sigma.
    Therefore K=sum_n f_n Tr{v[Gr sigma P_n+P_n sigma Ga]}.
    Avoid all interband Berry and geometry conventions in this calculation.
    Arrays must have trailing dimensions 2x2; works with arbitrary shape.
    """
    H=np.asarray(H,dtype=complex); velocity=np.asarray(velocity,dtype=complex)
    if H.shape[-2:]!=(2,2) or H.shape!=velocity.shape:raise ValueError('H,v must have common 2x2 shape')
    if not np.isfinite(omega) or omega<=0 or eta<0:raise ValueError('omega>0 and eta>=0 required')
    if spin not in SPINS:raise ValueError('Pauli axis invalid')
    e,U=np.linalg.eigh(H)
    f=(e<mu).astype(float)
    s=SPINS[spin]
    out=np.zeros(e.shape[:-1],dtype=complex)
    for n in range(2):
        u=U[..., :,n]; proj=u[..., :,None]*u.conj()[...,None,:]
        er=e[...,n]
        gr=np.linalg.inv((er+omega+1j*eta)[...,None,None]*EYE-H)
        ga=np.linalg.inv((er-omega-1j*eta)[...,None,None]*EYE-H)
        # Purely intraband terms from Eq.(44) cancel algebraically as
        # 1/(omega+i eta)+1/(-omega-i eta)=0. Remove *only* that
        # machine-cancellation, keeping the cross-band resolvent intact.
        q=EYE-proj
        response=q@gr@s@proj+proj@s@ga@q
        out+=f[...,n]*np.trace(velocity@response,axis1=-2,axis2=-1)
    return out

def wavepacket_cos_sin(H,velocity,spin,omega,mu=0.):
    """Per-k physical a^(1)(t)/b = Pcos*cos(omega*t)+Psin*sin(omega*t).

    Derived directly from SM Eq.(21),(24), with B(t)=B cos(omega t).
    Returns C=partial_k(<sigma>_occupied) from the *same* two-band
    eigenvectors via spectral identity, only for an insulator one lower band.
    Exact at sub-gap omega, no relaxation, no frequency approximation.
    """
    e,U=np.linalg.eigh(H)
    if not np.all((e[...,0]<mu)&(e[...,1]>mu)):raise ValueError('filled lower and empty upper required')
    gap=e[...,1]-e[...,0]
    if np.any(gap<=omega):raise ValueError('subgap drive required')
    v01=np.einsum('...i,...ij,...j->...',U[..., :,0].conj(),velocity,U[..., :,1])
    s10=np.einsum('...i,ij,...j->...',U[..., :,1].conj(),SPINS[spin],U[..., :,0])
    # Main Eq(5): r_01=i v_01/gap; A=r_01 sigma_10.
    A=1j*v01*s10/gap
    up=1/(omega-gap);um=1/(-omega-gap)
    # SM Eq(24): 2Re[M_10 r_01] = -b Re[A(e^-iwt up + e^+iwt um)].
    cos_amp=-np.real(A)*(up+um)
    sin_amp=-np.imag(A)*(up-um)
    # SM/main geometry: Z=-2 Im(A) and Q=Re(A).
    Z=-2*np.imag(A)
    return {'P_cos_per_b':cos_amp,'P_sin_per_b':sin_amp,
            'C_static':Z,'Q_main':np.real(A),'Z_main':Z,
            'K_from_wavepacket':Z-omega*sin_amp+1j*omega*cos_amp,
            'K_main_lowfrequency':Z+1j*omega*2*np.real(A)/gap,
            'P_sin_main_Eq10_per_b':omega*Z/gap**2,
            'P_sin_SM_Eq34_per_b':sin_amp} 

def analytical_cosine_harmonics(N=256):
    """Under integral dt exp(i wt), Bcos => pi*B*(delta+ + delta-).
    A normalized FFT of one exact cycle measures harmonic amplitude B/2.
    """
    samples=np.cos(2*np.pi*np.arange(N)/N)
    positive=np.fft.fft(samples)[1]/N
    negative=np.fft.fft(samples)[-1]/N
    return float(positive.real),float(negative.real)
