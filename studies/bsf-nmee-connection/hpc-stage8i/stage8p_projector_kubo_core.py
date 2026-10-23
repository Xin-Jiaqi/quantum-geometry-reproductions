"""Stage8p: one-convention source/Kubo-projector interface in a finite Bloch basis.

Zero-job data analysis; does not itself read WSe2 files. Exactly one common
spin-only source convention for all methods:
 H = diag(E) + lambda * (SPN_x/2) - mu_B * B * SPN_z
where SPN stores Pauli sigma components. No extra g=2 applied to sigma_z.

The resolvent and the projector finite-difference routes share a COVARIANT
k-projector derivative constructed from already-existing MMN and eigenframes
at both NEIGHBOR k points. This is not mathematically independent of those
MMN links and still suffers finite-96-band leakage; report it honestly.
"""
from __future__ import annotations
import numpy as np

MU_B=5.7883818060e-5


def h_source(eig,spn_x,spn_z,lam,B,magnetic_prefactor=MU_B):
    """spn_* are Pauli sigma (not sigma/2)."""
    E=np.asarray(eig,dtype=float)
    H=np.diag(E).astype(complex)+.5*float(lam)*spn_x
    H-=float(magnetic_prefactor)*float(B)*spn_z
    if not np.allclose(H,H.conj().T,rtol=0,atol=1e-8):
        raise ValueError("Non-Hermitian source Hamiltonian")
    return H


def frame_and_projector(H,occ):
    e,v=np.linalg.eigh(H)
    if not (0<occ<len(e)):raise ValueError("Need occupied and unoccupied states")
    gap=float(e[occ]-e[occ-1])
    if gap<=0:raise ValueError("Occupied direct gap is closed")
    u=v[:,:occ]
    return e,v,u@u.conj().T,gap


def analytic_source_derivative(e,v,spin_source,occ):
    """Analytic dP/dlambda from interband resolvent, source=Sx/2.

    dP_mn (m unocc,n occ) = S_mn/(e_n-e_m), plus h.c.
    """
    O=v.conj().T@spin_source@v
    A=np.zeros_like(O,dtype=complex)
    for n in range(occ):
        for m in range(occ,len(e)):
            den=e[n]-e[m]
            if abs(den)<1e-8:raise ValueError("Small occupied-empty denominator")
            A[m,n]=O[m,n]/den
            A[n,m]=O[n,m]/den
    return v@A@v.conj().T


def finite_source_derivative(eig,spn_x,spn_z,B,occ,dl,magnetic_prefactor=MU_B):
    if dl<=0:raise ValueError("Invalid spin-source step")
    def P(l):
        return frame_and_projector(h_source(eig,spn_x,spn_z,l,B,magnetic_prefactor),occ)[2]
    return (P(+dl)-P(-dl))/(2*dl)


def covariant_k_projector(eig_p,sx_p,sz_p,M_p,
                          eig_m,sx_m,sz_m,M_m,
                          B,occ,dk,magnetic_prefactor=MU_B):
    """Map neighbor projectors into the CENTRAL Bloch basis.

    M_p,m = <u(center)|u(k+/-)> in same 96-band finite basis.
    Pp_c = M_p Pp M_p^dagger; similarly Pm.
    No arbitrary neighbor band phases survive. In a truncated Bloch window,
    Pp_c may not be idempotent; report this leakage and link unitarity.
    """
    if dk<=0:raise ValueError("Invalid reciprocal Cartesian step")
    Hp=h_source(eig_p,sx_p,sz_p,0,B,magnetic_prefactor)
    Hm=h_source(eig_m,sx_m,sz_m,0,B,magnetic_prefactor)
    _,_,Pp,gp=frame_and_projector(Hp,occ)
    _,_,Pm,gm=frame_and_projector(Hm,occ)
    nb=Pp.shape[0]
    if M_p.shape!=(nb,nb) or M_m.shape!=(nb,nb):
        raise ValueError("MMN link dimensions do not match eig/spin window")
    if not np.all(np.isfinite(M_p)) or not np.all(np.isfinite(M_m)):
        raise ValueError("MMN links contain nonfinite values")
    pp=M_p@Pp@M_p.conj().T
    pm=M_m@Pm@M_m.conj().T
    I=np.eye(nb)
    leakage={
        "max_neighbour_projector_idempotency_error":
          float(max(np.linalg.norm(pp@pp-pp,ord="fro"),
                    np.linalg.norm(pm@pm-pm,ord="fro"))),
        "max_link_nonunitarity_fro":
          float(max(np.linalg.norm(M_p.conj().T@M_p-I,ord="fro"),
                    np.linalg.norm(M_m.conj().T@M_m-I,ord="fro"))),
        "min_neighbour_occ_empty_gap_eV":min(gp,gm)
    }
    D=(pp-pm)/(2*dk)
    if not np.allclose(D,D.conj().T,atol=1e-7):
        raise ValueError("Covariant k projector derivative is not Hermitian")
    return D,leakage


def mixed_geom(P,dPl,dPk):
    """i Tr P[dPl,dPk] with both derivatives Hermitian."""
    return float((1j*np.trace(P@(dPl@dPk-dPk@dPl))).real)


def compare_local_three_way(center,plus,minus,links,
                            *,B,occ,dk,dl,magnetic_prefactor=MU_B):
    """Inputs tuples (E,SPN_x,SPN_z) for center and +/- neighbor k.

    Returns resolvent-source vs projector-source F at fixed B; third-route
    Wilson must be compared outside this function with an independent
    plaquette engine. F is NOT yet a magnetic derivative Gamma.
    """
    e,sx,sz=center
    E,V,P,gap=frame_and_projector(h_source(e,sx,sz,0,B,magnetic_prefactor),occ)
    dPk,leak=covariant_k_projector(
        *plus,links[0],*minus,links[1],B,occ,dk,magnetic_prefactor)
    dPl_analytic=analytic_source_derivative(E,V,sx/2,occ)
    dPl_fd=finite_source_derivative(e,sx,sz,B,occ,dl,magnetic_prefactor)
    Fr=mixed_geom(P,dPl_analytic,dPk)
    Fp=mixed_geom(P,dPl_fd,dPk)
    return {
        "F_resolvent":Fr,"F_projector_fd":Fp,
        "absolute_source_derivative_discrepancy":abs(Fr-Fp),
        "occ_unocc_gap_eV":gap,
        "source_convention":"lambda*SPN_x/2",
        "magnetic_convention":"-mu_B*B*SPN_z (SPN_z=sigma_z)",
        "MMN_subspace_leakage":leak,
        "material_static_gamma_computed":False
    }
