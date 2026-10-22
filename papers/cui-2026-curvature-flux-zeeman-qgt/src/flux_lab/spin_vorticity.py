"""Band spin vorticity, Cui et al. 2026 main Eqs. (20)-(28).

Conventional spin O_a^c={dH/dk_a,sigma_c}/4, fixed k-independent
orthonormal basis and k-independent *physical* Pauli-spin operator.
Omega_bar^A_xy = sum_m (E_n-E_m) Omega^A_nm,xy = (curl O_n)_z/2.
All results are model natural units and require nondegenerate bands.
The spectral geometric RTA part is not generally the full periodic-lattice
RTA because explicit Hessian vertex can be nonzero.
"""
from __future__ import annotations
import numpy as np
from .spin_current import (_validated, spin_current_operator,
    spin_orbital_benchmark, spin_current_derivative)


def spin_vorticity_spectral(h, velocities, physical_spin):
    """Independent pair-weighted spin metric/vorticity and Eq.(28) reduced RTA.

    Returns shape (nbands,2,2) symmetric metric and antisymmetric matrix,
    plus vector z component, geometric RTA and intrinsic Eq.(128) response.
    """
    _,v,s,_,e,U=_validated(h,velocities,physical_spin)
    ob=np.einsum('in,aij,jm->anm',U.conj(),spin_current_operator(v,s),U,optimize=True)
    vb=np.einsum('in,aij,jm->anm',U.conj(),v,U,optimize=True)
    n=len(e)
    weighted=np.zeros((n,2,2),complex)
    unweighted=np.zeros((n,2,2),complex)
    for i in range(n):
        for j in range(n):
            if i==j:continue
            P=np.outer(ob[:,i,j],vb[:,j,i]);delta=e[i]-e[j]
            weighted[i]+=P/delta
            unweighted[i]+=P/delta**2
    gbar=((weighted+np.swapaxes(weighted,1,2))/2).real
    omega_bar=-(weighted-np.swapaxes(weighted,1,2)).real
    gA=((unweighted+np.swapaxes(unweighted,1,2))/2).imag
    omegaN=-(unweighted-np.swapaxes(unweighted,1,2)).imag
    intrinsic=-2*unweighted.imag
    reduced_RTA=-2*weighted.real
    return dict(energies=e,bar_gN=gbar,bar_omegaA=omega_bar,
                bar_omegaA_z=omega_bar[:,0,1],
                reduced_RTA=reduced_RTA,
                reduced_RTA_from_sectors=-2*gbar+omega_bar,
                gA=gA,omegaN=omegaN,intrinsic=intrinsic,
                intrinsic_from_sectors=-2*gA+omegaN)


def vorticity_with_hessian(h, velocities, hessian, physical_spin):
    """Eq.(27) plus full periodic band derivative; Hessian is symmetric in a,b.

    Antisymmetrizing the *full* derivative cancels the explicit Hessian.
    This cancellation does NOT justify omitting Hessian from symmetric RTA.
    """
    q=spin_vorticity_spectral(h,velocities,physical_spin)
    d=spin_current_derivative(h,velocities,hessian,physical_spin)
    # d['full_RTA'][n,a,b] = -partial_b <O_a>_n
    # (curl O)_z / 2 = [partial_x<O_y>-partial_y<O_x>]/2
    curl_half=(d['full_RTA'][:,0,1]-d['full_RTA'][:,1,0])/2
    vertex_curl=(d['hessian_vertex'][:,1,0]-d['hessian_vertex'][:,0,1])/2
    return {**q,'curl_half_from_full_derivative':curl_half,
            'vertex_curl_half':vertex_curl,
            'full_RTA':d['full_RTA'],
            'symmetric_hessian_correction':-d['hessian_vertex']}


def vorticity_finite_difference(kx,ky,*,periodic=True,step=1e-5):
    """Eigenstate-only 1/2(d_x<O_y>-d_y<O_x>), no QGT or Hessian used."""
    if step<=0:raise ValueError('finite-difference step must be positive')
    grads=np.zeros((2,4,2),float)
    for b in (0,1):
        values=[]
        for offset in (+step,-step):
            k=[kx,ky];k[b]+=offset
            h,v,_,s=spin_orbital_benchmark(*k,periodic=periodic)
            _,U=np.linalg.eigh(h)
            O=spin_current_operator(v,s)
            values.append(np.stack([np.diag(U.conj().T@oa@U).real for oa in O],axis=1))
        grads[b]=(values[0]-values[1])/(2*step)
    return (grads[0,:,1]-grads[1,:,0])/2
