"""Stage7i: gauge-invariant OCCUPIED-SUBSPACE Wilson plaquettes for static NMEE.

Use a rectangular Wilson loop in parameter space (Cartesian k_j, lambda*S_i)
at each k, then take a finite Bz derivative. All four loop vertices have
the entire occupied subspace; arbitrary U(Nocc) rotations cancel around it.

This is a validation of the gapped zero-temperature EQUILIBRIUM static
mixed Berry curvature, not the paper's nonequilibrium driven-DC NMEE.
For a nontrivial Wannier basis connection the link matrices must include
basis overlaps. This tight-binding module uses the original absolute-orbital
Bloch embedding; DO NOT feed hr.dat matrices alone into its overlap routine.
"""
from dataclasses import replace
import numpy as np
from nmee_lab.clean_core.core import MU_B
from .sliding_tight_binding import SlidingBilayer, SlidingSpec


def _occupied(model, k, source_matrix, field_T, source_eV):
    H=model.matrices(k)
    hh=H + source_eV*source_matrix - (
        model.mono.g_factor*MU_B*field_T*model.spins["total_z"])
    eigenvalues, eigenvectors=np.linalg.eigh(hh)
    mu=model.spec.mu_eV
    dist=min(float(np.min(eigenvalues[:,4]-mu)),
             float(np.min(mu-eigenvalues[:,3])))
    if dist <= 0:
        raise ValueError("Source plaquette closes the occupied/empty gap")
    return eigenvectors[:,:,:4],dist


def _field_loop_phase(model,k,spin_axis,k_axis,field_T,dk,dl):
    if dk<=0 or dl<=0:
        raise ValueError("Positive parameter loop steps required")
    dr=np.eye(2)[k_axis]*dk/2
    spin=model.spins["total_"+spin_axis]
    corners=((k-dr,-dl/2),(k+dr,-dl/2),
             (k+dr,+dl/2),(k-dr,+dl/2))
    frames=[]; gap=float("inf")
    for kp,lam in corners:
        frame,m=_occupied(model,kp,spin,field_T,lam)
        frames.append(frame);gap=min(gap,m)
    W=np.ones(len(k),complex)
    for a,b in zip(frames,frames[1:]+frames[:1]):
        O=np.einsum("nai,naj->nij",a.conj(),b,optimize=True)
        det=np.linalg.det(O)
        W *= det/(abs(det)+1e-300)
    return np.angle(W),gap


def wilson_static_kernel(model, k, *, spin_axis="x", k_axis=0,
                         step_T=.25,dk=.005,dl=.000125):
    """Local d_B F_(lambda_i,k_j) for an Nx2 k array.

    lambda multiplies Pauli S_i in eV. The Wilson oriented loop is
    k+,lambda+,k-,lambda-.  Its phase/(dk*dl) equals F_lambda,k.
    """
    if spin_axis not in ("x","y") or k_axis not in (0,1):
        raise ValueError("Only in-plane tensor indices are supported")
    if step_T<=0:raise ValueError("positive Zeeman field step required")
    k=np.atleast_2d(np.asarray(k,float))
    if k.shape[1:]!=(2,):raise ValueError("k is (N,2) Cartesian nm^-1")
    plus,gp=_field_loop_phase(model,k,spin_axis,k_axis,+step_T,dk,dl)
    minus,gm=_field_loop_phase(model,k,spin_axis,k_axis,-step_T,dk,dl)
    return (np.angle(np.exp(1j*(plus-minus)))/(dk*dl*2*step_T),
            min(gp,gm))


def adaptive_valley_disks(model, *, radius=.5, radial_n=22, angular_n=60):
    if radius<=0 or radial_n<2 or angular_n<4:
        raise ValueError("Invalid disk integration")
    reciprocal=model.model().reciprocal
    centers=np.array([[1/3,2/3],[2/3,1/3]])@reciprocal
    roots,gw=np.polynomial.legendre.leggauss(radial_n)
    rr=radius*(roots+1)/2
    rw=radius*gw/2
    phi=2*np.pi*(np.arange(angular_n)+.5)/angular_n
    q=(rr[:,None,None]*
       np.stack((np.cos(phi),np.sin(phi)),axis=-1)[None,:,:]).reshape(-1,2)
    w=np.repeat(rr*rw,angular_n)/(2*np.pi*angular_n)
    return np.concatenate((q+centers[0],q+centers[1])),np.tile(w,2)


def compare_valley_wilson_kubo(spec=SlidingSpec(), *, spin_axis="x",
                              k_axis=0, radius=.5,radial_n=22,angular_n=60,
                              dk=.005,dl=.000125,step_T=.25,batch=512):
    """Compare independent Wilson occupied-overlap vs static spectral kernel.

    NO full-BZ claim. Adaptive disks only. Physical spin and velocity are
    evaluated from the same original microscopic Hamiltonian.
    """
    from .valley_kp import static_kernel
    model=SlidingBilayer(spec)
    k,w=adaptive_valley_disks(model,radius=radius,radial_n=radial_n,
                              angular_n=angular_n)
    wtotal=0.;ktotal=0.;mingap=float("inf")
    idx=0 if spin_axis=="x" else 1
    for start in range(0,len(k),batch):
        kk=k[start:start+batch]
        ww=w[start:start+batch]
        kernel,gap=wilson_static_kernel(
            model,kk,spin_axis=spin_axis,k_axis=k_axis,
            step_T=step_T,dk=dk,dl=dl)
        H=model.matrices(kk)
        V=np.stack([model.matrices(kk,axis=j) for j in range(2)])
        kb=static_kernel(H,V,[model.spins["total_x"],model.spins["total_y"]],
                         model.spins["total_z"],spec.mu_eV,
                         g_factor=model.mono.g_factor,step_T=step_T)
        wtotal+=float(ww@kernel)
        ktotal+=float(ww@kb[:,idx,k_axis])
        mingap=min(mingap,gap)
    return {"wilson":wtotal,"kubo":ktotal,"difference":wtotal-ktotal,
        "relative_difference":abs(wtotal-ktotal)/max(abs(ktotal),1.e-100),
        "minimum_sampled_mu_gap_eV":mingap,
        "u_fractional":list(spec.shift),"tunneling_eV":spec.tunneling_eV,
        "spin_axis":spin_axis,"k_axis":k_axis,"dk_nm_inverse":dk,
        "lambda_source_eV":dl,"B_step_T":step_T,"disk_radius_nm_inverse":radius,
        "radial_n":radial_n,"angular_n":angular_n,
        "scope":"T=0 insulating STATIC response on two valley disks, not driven DC"}
