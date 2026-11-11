"""Stage 7d: direct valley chirality and Gaussian hopping/velocity identity.

Exact within the declared orbital embedding and finite Gaussian hopping shell.
The resulting chiral hopping is NOT by itself a sign theorem for NMEE.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.models import SZ
from .sliding_tight_binding import SlidingBilayer

SPIN_PI_Z = np.kron(np.eye(4, dtype=complex), SZ)


def orbital_hopping_vertices(bilayer: SlidingBilayer, k):
    """Return 2x2 hopping T, momentum V_j, and physical sliding J_j.

    For D=R+r_b-r_a+d, t_D~exp(-D²/xi²) and Bloch exp(i k.D):
      dT/dd_j = i k_j T + (2i/xi²) dT/dk_j.
    The equality is gauge- and hopping-law-dependent (not universal).
    """
    k=np.asarray(k,float)
    if k.shape != (2,):
        raise ValueError("k must be 2D Cartesian reciprocal nm")
    xi=bilayer.spec.decay_nm
    T=np.zeros((2,2),complex)
    V=np.zeros((2,2,2),complex)
    J=np.zeros((2,2,2),complex)
    for ia,ib,displacements,amplitudes in bilayer.hops:
        weights=np.exp(1j*displacements@k)*amplitudes
        T[ia,ib]+=np.sum(weights)
        for j in range(2):
            V[j,ia,ib]+=np.sum(1j*displacements[:,j]*weights)
            J[j,ia,ib]+=np.sum((1j*k[j]-2*displacements[:,j]/xi**2)*weights)
    return T,V,J


def chiral_valley_summary(spec):
    """Check at K/K' for the actual 8-band source, not fitted k.p data."""
    obj=SlidingBilayer(spec)
    g=obj.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@g
    rows=[]
    for nu,center in enumerate(centers):
        T,V,J=orbital_hopping_vertices(obj,center)
        source=obj.tunneling(center[None,:])[0]
        source_vel=[obj.tunneling(center[None,:],axis=j)[0] for j in range(2)]
        err=max(np.max(np.abs(T-source[::2,::2])),
                max(np.max(np.abs(V[j]-source_vel[j][::2,::2])) for j in range(2)))
        identity=np.max(np.abs(J-1j*center[:,None,None]*T[None,:,:]
                               -2j/spec.decay_nm**2*V))
        rows.append(dict(valley=nu,center=center.tolist(),
                         hopping_AB_meV=float(abs(T[0,1])*1e3),
                         hopping_AA_meV=float(abs(T[0,0])*1e3),
                         velocity_AB_x_meVnm=float(abs(V[0,0,1])*1e3),
                         vy_over_vx_AB=[float(x) for x in
                                        (np.real(V[1,0,1]/V[0,0,1]),
                                         np.imag(V[1,0,1]/V[0,0,1]))],
                         vy_over_vx_BA=[float(x) for x in
                                        (np.real(V[1,1,0]/V[0,1,0]),
                                         np.imag(V[1,1,0]/V[0,1,0]))],
                         microscopic_velocity_residual=float(err),
                         slide_velocity_identity_residual=float(identity)))
    return rows


def rashba_spin_rotation_check(spec,k):
    """U_z H(+alpha) U_z^dagger=H(-alpha), Sx odd, Sz even."""
    pos=SlidingBilayer(spec)
    neg=SlidingBilayer(replace(spec,rashba_soc_eV=-spec.rashba_soc_eV))
    x=np.asarray(k,float)[None,:]
    U=SPIN_PI_Z
    H=pos.matrices(x)[0]
    Hneg=neg.matrices(x)[0]
    v=pos.matrices(x,axis=0)[0]
    vneg=neg.matrices(x,axis=0)[0]
    Sx=pos.spins["total_x"]
    Sz=pos.spins["total_z"]
    return {
        "Hamiltonian_residual":float(np.max(np.abs(Hneg-U@H@U.conj().T))),
        "velocity_residual":float(np.max(np.abs(vneg-U@v@U.conj().T))),
        "Sx_odd_residual":float(np.max(np.abs(U@Sx@U.conj().T+Sx))),
        "Sz_even_residual":float(np.max(np.abs(U@Sz@U.conj().T-Sz))),
    }
