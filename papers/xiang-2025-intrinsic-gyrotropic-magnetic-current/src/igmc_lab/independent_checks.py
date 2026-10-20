"""Independent two-band IGMC benchmarks, not reusing Zeeman QGT band sums.

1. Finite-field Wilson-loop Berry phase (insulating polarization derivative).
2. Direct orbital-basis Liouville matrix inversion at finite frequency.

All k are 1/Angstrom, H/b/omega energy in eV, hbar=e=1 here.
"""
from __future__ import annotations
import numpy as np
from .core import dirac_grid,dirac_hamiltonian,SPINS


def finite_field_wilson_derivative(n: int, *, b_eV: float=1e-4, a:float=1.,
                                   beta:float=.25, v:float=.825, M:float=.1,
                                   spin: str='x') -> float:
    """d/d(b) [int_BZ A_x(k) d^2k/(2pi)^2], H(b)=H(0)-b*sigma_spin.

    Use occupied-band Wilson-loop link products only (NO velocity matrix,
    Berry interband connection, Zeeman kernel, or analytic band sum).
    The BZ torus phase is computed for every fixed k_y. Finite differences
    of principal branch phases are valid if each change stays below pi.
    """
    if n<9 or b_eV<=0 or not np.isfinite(b_eV): raise ValueError('invalid grid or magnetic energy')
    kx,ky=dirac_grid(n,a)
    H,_=dirac_hamiltonian(kx,ky,beta=beta,a=a,v=v,M=M)
    op=SPINS[spin]
    phi=[]
    for bb in (-b_eV,b_eV):
        eig,states=np.linalg.eigh(H-bb*op)
        if np.any(eig[...,1]-eig[...,0]<1e-8):raise ValueError('gap closure at finite field')
        occupied=states[..., :,0] # (kx,ky,orbital)
        adjacent=np.roll(occupied,-1,axis=0)
        overlap=np.sum(occupied.conj()*adjacent,axis=-1)
        if np.any(abs(overlap)<1e-9):raise ValueError('Wilson link ill conditioned')
        phi.append(-np.angle(np.prod(overlap/abs(overlap),axis=0)))
    # Circular phase difference modulo 2pi is immune to Wilson winding gauge.
    dphi=np.angle(np.exp(1j*(phi[1]-phi[0])))
    # dk_y=2pi/(Na). Wilson phase is the kx integral; divide by (2pi)^2.
    return float(np.sum(dphi)/(2*b_eV*n*2*np.pi*a))


def orbital_liouville_kernel(n:int, *, mu:float=0., omega_eV:float=.004135667694,
                              a:float=1., beta:float=.25,v:float=.825,M:float=.1,
                              alpha:str='x', spin:str='x') -> complex:
    """Compute K(omega)=-Tr[v_x delta_rho]/b, via orbital Liouvillian solve.

    Solve [omega I - (H. - .H)] rho^(1) = -[sigma, rho0], where the
    magnetic perturbation is -b*sigma exp(-i*omega*t). No eigenvector
    interband formula in this response or geometric normalization.
    Requires zero-temperature insulator with chemical potential in gap.
    """
    if n<9:raise ValueError('grid too small')
    kx,ky=dirac_grid(n,a)
    H,dH=dirac_hamiltonian(kx,ky,a=a,beta=beta,v=v,M=M)
    energy=np.sqrt(np.real(np.trace(H@H,axis1=-2,axis2=-1))/2)
    if np.any(energy<=abs(mu)) or omega_eV>=float(2*energy.min()):
        raise ValueError('only gapped nonresonant T=0 occupations supported')
    I=np.eye(2,dtype=complex)
    rho0=(I-H/energy[...,None,None])/2
    s=SPINS[spin]
    rhs=-(s@rho0-rho0@s)
    # Vectorization in C order (ij) => H_ik delta_jl - delta_ik H_lj.
    L=np.empty(H.shape[:-2]+(4,4),complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    L[...,2*i+j,2*k+l] = (
                        omega_eV*(i==k)*(j==l)
                        -H[...,i,k]*(j==l) + (i==k)*H[...,l,j]
                    )
    d_rho=np.linalg.solve(L,rhs.reshape(H.shape[:-2]+(4,1)))[...,0].reshape(H.shape)
    current_per_b=np.trace(dH[alpha]@d_rho,axis1=-2,axis2=-1)
    # Eq.(8) geometric K = - Tr(v rho^(1))/b.
    return complex(-np.mean(current_per_b)/a**2)
