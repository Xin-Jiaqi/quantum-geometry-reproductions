"""Independent eigenvector-free two-band Zeeman current identities.

Derive both current kernels from the normalized Bloch-sphere vector
n=d/|d| and Pauli trace identities. Does not import the published
spectral/QGT/Kubo implementation. Frequencies here are photon energies,
hbar*angular_frequency, in eV; k vectors have inverse-Angstrom units.

H(k)=d_0(k) I + d(k).sigma, Zeeman drive H1=-b sigma_beta.
For ordered lower/upper bands (gap=2E):
  Z^alpha_beta(lower,upper) = -partial_kalpha n_beta
  [2 Q^alpha_beta / gap](lower,upper)
       = d.dot(partial_kalpha d cross e_beta) / (2 E^3)
The definitions equal the main paper Eqs. (5), (6), (11), (12).
The finite-frequency pair kernel for nonresonant real omega is
  K = (f_lower-f_upper) (Z+i*omega*Dperomega)/(1-(omega/gap)^2).
The two quadratures have different TR/domain parity. This is NOT an
ad-hoc factor-two repair of the conflicting SM Eq.(34) intermediate.
"""
from __future__ import annotations
import numpy as np

_E = {"x":np.array((1.,0.,0.)), "y":np.array((0.,1.,0.)), "z":np.array((0.,0.,1.))}


def dvector_dirac(kx,ky,*, beta=.25, a=1.,v=.825,M=.1,mass_sign=1):
    """Return d, (d_x', d_y') for the PRL Eq.(17) lattice Dirac model.

    Signs and lattice factors are transcribed directly from the paper;
    mass_sign=-1 denotes the time-reversed magnetic domain.
    """
    if a<=0 or mass_sign not in (-1,1):raise ValueError("invalid lattice or mass sign")
    x,y=a*np.asarray(kx),a*np.asarray(ky)
    base=np.broadcast_arrays(x,y)
    x,y=base
    zero=np.zeros_like(x)
    d=np.stack((v/a*np.sin(y),-v/a*np.sin(x),
               mass_sign*(beta+2*M/a**2*(2-np.cos(x)-np.cos(y)))),axis=-1)
    grad={
        'x':np.stack((zero,-v*np.cos(x),mass_sign*2*M/a*np.sin(x)),axis=-1),
        'y':np.stack((v*np.cos(y),zero,mass_sign*2*M/a*np.sin(y)),axis=-1),
    }
    return d,grad


def pauli_trace_geometry(d,grad,alpha='x',beta='x',gap_tol=1.e-10):
    """Return C (Zeeman curvature), D/omega, and gap without diagonalizing H.

    Each result is per lower->upper pair and is well-defined only if E>0.
    """
    if beta not in _E or alpha not in grad:raise ValueError('unsupported index')
    d=np.asarray(d,dtype=float)
    v=np.asarray(grad[alpha],dtype=float)
    if d.shape!=v.shape or d.shape[-1]!=3:raise ValueError('expected matching Cartesian d and derivative')
    E=np.linalg.norm(d,axis=-1)
    if np.any(~np.isfinite(E)) or np.any(2*E<=gap_tol):
        raise ValueError('nondegenerate finite-energy model required')
    e=_E[beta]
    dn=np.sum(d*v,axis=-1)
    # Z=-partial_{k_alpha}(d_beta/E) = -derivative of the upper-band spin.
    Z=-v[...,"xyz".index(beta)]/E + d[...,"xyz".index(beta)]*dn/E**3
    D=np.sum(d*np.cross(v,e),axis=-1)/(2*E**3)
    return {'Z01':Z,'D_per_omega01':D,'gap':2*E}


def occupation_zeroT_from_energy(E,mu):
    """Occupation difference lower minus upper for +/-E bands."""
    return (np.asarray(mu)>-E).astype(float)-(np.asarray(mu)>E).astype(float)


def aggregate_response(geometric,mu=0.,photon_eV=.004135667694,*,a=1.):
    """BZ midpoint quadrature for a fully specified two-band dispersion."""
    if a<=0:raise ValueError('positive lattice length required')
    E=geometric['gap']/2
    if photon_eV<0 or photon_eV>=float(np.min(geometric['gap'])):
        raise ValueError('requires real nonresonant photon energy')
    delta=occupation_zeroT_from_energy(E,mu)
    den=1-(photon_eV/geometric['gap'])**2
    C=np.mean(delta*geometric['Z01'])/a**2
    D=np.mean(delta*geometric['D_per_omega01'])/a**2
    # Independent closed-form matrix-algebra Kubo, with full omega denominators.
    K=np.mean(delta*(geometric['Z01']+1j*photon_eV*geometric['D_per_omega01'])/den)/a**2
    return {'sigma_C':float(C),'sigma_D_per_omega':float(D),
            'sigma_D':float(D*photon_eV),
            'Kubo_real':float(K.real),'Kubo_imag':float(K.imag)}


def midpoint_dirac(N,mu=0.,photon_eV=.004135667694,*,alpha='x',spin='x',**kwargs):
    """Compute the paper's Dirac response directly from 3-vectors, no eigenstates."""
    if N<3:raise ValueError('N >=3')
    a=kwargs.get('a',1.)
    k=(-np.pi+(np.arange(N)+.5)*2*np.pi/N)/a
    x,y=np.meshgrid(k,k,indexing='ij')
    d,der=dvector_dirac(x,y,**kwargs)
    z=pauli_trace_geometry(d,der,alpha,spin)
    return aggregate_response(z,mu,photon_eV,a=a)
