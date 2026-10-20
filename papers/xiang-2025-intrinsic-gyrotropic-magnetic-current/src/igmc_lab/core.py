"""Paper-convention Zeeman QGT and frequency-resolved IGMC for exactly two nondegenerate bands.

All H in eV, k in 1/Angstrom, dH/dk in eV Angstrom. hbar=e=1 for
algebraic response coefficients. Physical SI current conversion unverified.
Degenerate momenta raise, never broaden by default.
"""
from __future__ import annotations
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0,1],[1,0]],complex)
SY = np.array([[0,-1j],[1j,0]],complex)
SZ = np.array([[1,0],[0,-1]],complex)
SPINS={'x':SX,'y':SY,'z':SZ}
KB_EV_K=8.617333262145e-5
# 2022 CODATA (via scipy constants), frozen in source for no-dependency offline execution.
SI_E_C=1.602176634e-19
SI_HBAR_J_S=1.054571817e-34
SI_MU_B_J_T=9.2740100783e-24
SI_ANGSTROM_M=1.e-10

def energy_from_frequency(value, *, protocol):
    """Compute photon energy hbar*omega in eV with an EXPLICIT frequency convention.

    protocol='cyclic_hz': value=f in cycles/s, omega=2*pi*f.
    protocol='angular_rad_s': value=omega in radians/s.
    Passing just 'Hz' is intentionally forbidden to avoid 2*pi ambiguity.
    """
    if not np.isfinite(value) or value < 0:
        raise ValueError('frequency must be a finite nonnegative number')
    if protocol == 'cyclic_hz':
        omega=value*2*np.pi
    elif protocol == 'angular_rad_s':
        omega=value
    else:
        raise ValueError('frequency protocol must be cyclic_hz or angular_rad_s')
    return float(SI_HBAR_J_S*omega/SI_E_C)

def physical_response_si(model_kernel_inv_angstrom, *, g_factor=2.0):
    """Current per magnetic field in A / (m T) from spin-geometric kernel [Angstrom^-1].

    Returns (-g*mu_B*sigma_C/D) in the published plotting convention,
    with electron charge and hbar restored. The minus is in the paper's
    definition of J, and the positive result here is the plotted -gmu_B*sigma.
    No empirical rescaling. g must be documented as a chosen input.
    """
    if not np.isfinite(g_factor) or g_factor < 0:
        raise ValueError('g_factor must be finite and nonnegative')
    return float(model_kernel_inv_angstrom *
                 (SI_E_C*g_factor*SI_MU_B_J_T/SI_HBAR_J_S) / SI_ANGSTROM_M)


def _rotate(mat, U):
    return np.matmul(np.matmul(U.conj().swapaxes(-2,-1),mat),U)

def occupation(eps, mu, temp_k):
    if temp_k < 0: raise ValueError('T cannot be negative')
    if temp_k == 0:
        return np.asarray(eps < mu,dtype=float)
    x=np.clip((eps-mu)/(KB_EV_K*temp_k),-700,700)
    return 1/(1+np.exp(x))

def _h_from_d(energy0, dx,dy,dz):
    h=np.empty(np.shape(dx)+(2,2), complex)
    h[...,0,0]=energy0+dz
    h[...,1,1]=energy0-dz
    h[...,0,1]=dx-1j*dy
    h[...,1,0]=dx+1j*dy
    return h

def dirac_grid(N, a=1.0):
    if N<3: raise ValueError('N must be >=3')
    k=(-np.pi+(np.arange(N)+0.5)*(2*np.pi/N))/a
    return np.meshgrid(k,k,indexing='ij')

def dirac_hamiltonian(kx,ky,*,beta=.25,a=1.,v=.825,M=.1,mass_sign=1):
    """Eq.(17) and analytic derivatives in eV and eV.Angstrom."""
    if a<=0: raise ValueError('a>0')
    if mass_sign not in (1,-1):
        raise ValueError('mass_sign=+1 is the printed Dirac model, -1 is its TR partner')
    ax,ay=a*kx,a*ky
    dx=v/a*np.sin(ay)
    dy=-v/a*np.sin(ax)
    dz=mass_sign*(beta+2*M/a**2*(2-np.cos(ax)-np.cos(ay)))
    H=_h_from_d(0,dx,dy,dz)
    dHx=_h_from_d(0,np.zeros_like(dx),-v*np.cos(ax),mass_sign*2*M/a*np.sin(ax))
    dHy=_h_from_d(0,v*np.cos(ay),np.zeros_like(dx),mass_sign*2*M/a*np.sin(ay))
    return H,{'x':dHx,'y':dHy}

def rashba_hamiltonian(kx,ky,*,lam=.6,mass_ratio=1.1,hbar2_over_2me=3.80998212):
    alpha=hbar2_over_2me/mass_ratio
    k2=kx*kx+ky*ky
    H=_h_from_d(alpha*k2,lam*ky,-lam*kx,np.zeros_like(kx))
    dHx=_h_from_d(2*alpha*kx,np.zeros_like(kx),-lam*np.ones_like(kx),np.zeros_like(kx))
    dHy=_h_from_d(2*alpha*ky,lam*np.ones_like(kx),np.zeros_like(kx),np.zeros_like(kx))
    return H,{'x':dHx,'y':dHy}

def pair_kernels(H,dH,alpha,beta,*,gap_tol=1e-10):
    """Return per-k pair geometric kernels, diagonal velocities/spin & amplitudes.

    z_01 = r_01^alpha sigma_10^beta with r_01=i v_01/(e1-e0).
    Distinct-band bands ordered by np.linalg.eigh.
    """
    if alpha not in dH or beta not in SPINS:raise ValueError('unsupported Cartesian direction')
    e,U=np.linalg.eigh(H)
    delta=e[...,1]-e[...,0]
    if np.any(delta<=gap_tol):raise ValueError('degenerate pair: nondegenerate Zeeman geometry undefined')
    V=_rotate(dH[alpha],U)
    S=_rotate(SPINS[beta],U)
    A=V[...,0,1]*S[...,1,0]
    z=1j*A/delta
    return {'energies':e,'gap':delta,'Z01':-2*np.imag(z),'Q01':np.real(z),
            'v_diag':V.diagonal(axis1=-2,axis2=-1).real,
            's_diag':S.diagonal(axis1=-2,axis2=-1).real,
            'amplitude01':A}

def kernel_response(pair,mu,temperature_K=0.,omega_eV=.001,*,weight=1.):
    if omega_eV<0:raise ValueError('nonnegative photon energy required')
    e=pair['energies'];gap=pair['gap'];A=pair['amplitude01'];Q=pair['Q01'];Z=pair['Z01']
    f=occupation(e,mu,temperature_K)
    df=f[...,0]-f[...,1]
    c=weight*np.mean(df*Z)
    dperomega=weight*np.mean(df*2*Q/gap)
    if omega_eV>=np.min(gap):
        raise ValueError('omega >= direct gap: nonresonant Kubo approximation not permitted')
    # Eq.(8), paired n=0,1 terms, eta=0; no phenomenological broadening.
    kubo=weight*np.mean(df*(A/(omega_eV-gap)-np.conj(A)/(omega_eV+gap)))
    return {'sigma_C':float(c),'sigma_D_per_omega':float(dperomega),
            'sigma_D':float(dperomega*omega_eV),
            'Kubo_real':float(kubo.real),'Kubo_imag':float(kubo.imag),
            'direct_gap_min_eV':float(np.min(gap))}

def dirac_kernel(kx,ky,alpha='x',beta='y',**model):
    H,dH=dirac_hamiltonian(kx,ky,**model)
    return pair_kernels(H,dH,alpha,beta)

def rashba_kernel(kx,ky,alpha='x',beta='y',**model):
    H,dH=rashba_hamiltonian(kx,ky,**model)
    return pair_kernels(H,dH,alpha,beta)

def dirac_response(N,mu=0.,omega_eV=.001,temperature_K=0.,alpha='x',beta='y',**model):
    a=model.get('a',1.)
    kx,ky=dirac_grid(N,a)
    pair=dirac_kernel(kx,ky,alpha,beta,**model)
    return kernel_response(pair,mu,temperature_K,omega_eV,weight=1/a**2)

def fermi_surface_response(pair,mu,temperature_K,*,weight=1.):
    """Independent Eq.(13) Fermi-surface route; requires T>0 for smooth f' kernel."""
    if temperature_K<=0:raise ValueError('finite temperature required for a smooth delta kernel')
    f=occupation(pair['energies'],mu,temperature_K)
    fprime=-f*(1-f)/(KB_EV_K*temperature_K)
    return float(-weight*np.mean(np.sum(fprime*pair['v_diag']*pair['s_diag'],axis=-1)))

def rashba_analytic(mu,*,lam=.6,mass_ratio=1.1,hbar2_over_2me=3.80998212):
    """Eq. (16) after removing the paper's overall e/hbar; units inverse Angstrom."""
    alpha=hbar2_over_2me/mass_ratio
    minimum=-lam**2/(4*alpha)
    if mu<minimum: return 0.0
    if mu<0: return float(np.sqrt(lam**2+4*alpha*mu)/(4*np.pi*alpha))
    return float(lam/(4*np.pi*alpha))

def rashba_integrate(mu,*,nrad=1200,nphi=96,kmax=.5,lam=.6,mass_ratio=1.1):
    """Polar BZ-free quadrature for Eq.(11), angular/radial midpoint grid.

    Need kmax above all Fermi surfaces; never use its value as a fit parameter.
    """
    radii=(np.arange(nrad)+.5)*(kmax/nrad)
    phi=(np.arange(nphi)+.5)*2*np.pi/nphi
    kx=radii[:,None]*np.cos(phi)[None,:]
    ky=radii[:,None]*np.sin(phi)[None,:]
    pair=rashba_kernel(kx,ky,lam=lam,mass_ratio=mass_ratio)
    f=occupation(pair['energies'],mu,0.)
    val=(f[...,0]-f[...,1])*pair['Z01']
    integral=np.sum(val*radii[:,None])*(kmax/nrad)*(2*np.pi/nphi)/(2*np.pi)**2
    return float(integral)
