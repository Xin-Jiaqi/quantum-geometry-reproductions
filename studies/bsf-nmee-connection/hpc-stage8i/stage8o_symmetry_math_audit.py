"""Stage8o independent mathematical audit for Stage8n symmetry & Wilson scaling.

Offline only; no HPC, no user raw files. All assertions about actual VASP
source data must remain qualified until operator sewing is tested on site.
"""
from __future__ import annotations
import math
import numpy as np

A_HEX = 3.319006624468
N_MESH = 21
OLD_CART_MIRROR = np.array([[0., -1.], [-1., 0.]])
VALID_FRAC_MIRROR = np.array([[0, 1], [1, 0]], dtype=int)


def direct_basis(a=A_HEX):
    """Columns are a1=(a,0) and a2=(-a/2,+sqrt(3)*a/2)."""
    return np.array([[a, -a/2.], [0., math.sqrt(3.)*a/2.]])


def reciprocal_basis(a=A_HEX):
    return 2.0 * math.pi * np.linalg.inv(direct_basis(a)).T


def fractional_action_of_cartesian(cart_R, a=A_HEX):
    A=direct_basis(a)
    return np.linalg.solve(A, np.asarray(cart_R,dtype=float)@A)


def is_crystal_lattice_rotation(cart_R, a=A_HEX, tol=1e-9):
    frac=fractional_action_of_cartesian(cart_R,a)
    U=np.rint(frac).astype(int)
    return {
        "valid":bool(np.max(np.abs(frac-U))<tol and
            abs(round(np.linalg.det(U)))==1),
        "fractional_direct_matrix":frac.tolist(),
        "nearest_integer_matrix":U.tolist(),
        "max_noninteger_residual":float(np.max(np.abs(frac-U)))
    }


def cartesian_action_of_fractional_direct(U,a=A_HEX):
    A=direct_basis(a)
    return A@np.asarray(U,dtype=float)@np.linalg.inv(A)


def fractional_k_map(u,v,U):
    """Crystallographic direct fractional r' = U r => k' = U^-T k."""
    U=np.asarray(U,dtype=int)
    Ui=np.linalg.inv(U).T
    q=Ui@np.array([u,v],float)
    return tuple(float(x%1) for x in q)


def discrete_mesh_permutation(U,n=N_MESH):
    """Bijection of all n² uniform reciprocal fractional points."""
    U=np.asarray(U,dtype=int)
    if abs(round(np.linalg.det(U)))!=1:
        raise ValueError("not unimodular crystallographic operation")
    R=np.rint(np.linalg.inv(U).T).astype(int)
    if not np.allclose(R,np.linalg.inv(U).T,atol=1e-12):
        raise ValueError("non-integer reciprocal map")
    pairs={}
    for i in range(n):
        for j in range(n):
            p=R@np.array([i,j],dtype=int)
            pairs[(i,j)]=(int(p[0]%n),int(p[1]%n))
    if len(set(pairs.values()))!=n*n:
        raise AssertionError("fractional map not a bijection")
    return pairs


def axial_spin_rotation(R2):
    """Axial 3D spin under vertical in-plane mirror R2; z coordinate fixed."""
    R3=np.eye(3)
    R3[:2,:2]=np.asarray(R2,dtype=float)
    return np.linalg.det(R3)*R3


def spin_gram(spins):
    """Real Hermitian spin-component Gram: C_ab=Re Tr[S_a^dagger S_b]."""
    if len(spins)!=3:
        raise ValueError("need x/y/z spin matrices")
    return np.array([[float(np.trace(a.conj().T@b).real) for b in spins]
                     for a in spins])


def rotated_spin_gram(C,R2):
    Q=axial_spin_rotation(R2)
    return Q@np.asarray(C)@Q.T


def reported_curvature_to_raw_phase(curvature,dk,dl):
    """Stage8m plaquette_phase already computes angle(det)/dk/dl."""
    return float(curvature*dk*dl)


def diagnose_stage8n_phase_example(a=A_HEX, n=N_MESH, dl=1e-4):
    # This sample is reproduced from Stage8n AB diag4 sanitized ledger.
    f_plus=0.004314914511363606
    f_minus=0.004315389651252636
    dk=float(np.linalg.norm(reciprocal_basis(a)[:,0])/n)
    angle_plus=reported_curvature_to_raw_phase(f_plus,dk,dl)
    angle_minus=reported_curvature_to_raw_phase(f_minus,dk,dl)
    return {
        "dk_Ang_inv":dk,
        "lambda_step_eV":dl,
        "loop_phase_dimensionless_plus":angle_plus,
        "loop_phase_dimensionless_minus":angle_minus,
        "raw_loop_phase_difference_abs":abs(angle_plus-angle_minus),
        "raw_difference_over_machine_epsilon":abs(angle_plus-angle_minus)/np.finfo(float).eps,
        "original_report_uses_divided_curvature_as_undivided_phase":True,
        "true_precision_claim":"Recompute numerical condition including products of 36x36 polar-unitary links; epsilon ratio alone is insufficient"
    }


def audit_summary():
    invalid=is_crystal_lattice_rotation(OLD_CART_MIRROR)
    valid_cart=cartesian_action_of_fractional_direct(VALID_FRAC_MIRROR)
    valid=is_crystal_lattice_rotation(valid_cart)
    return {
        "stage":"8o","source_stage":"8n",
        "old_reported_cartesian_mirror":OLD_CART_MIRROR.tolist(),
        "old_cartesian_mirror_preserves_actual_hex_lattice":invalid,
        "example_actual_fractional_crystal_mirror":VALID_FRAC_MIRROR.tolist(),
        "example_cartesian_crystal_mirror":valid_cart.tolist(),
        "example_is_valid_hex_crystal_mirror":valid,
        "full_21x21_mesh_symmetry_mapping_pairs":len(discrete_mesh_permutation(VALID_FRAC_MIRROR)),
        "horizontal_glide_inplane_k_map":"identity for pure m_z with in-plane t; nonsymmorphic sewing phase still required",
        "axial_spin_rule":"Q=det(R3) R3; compare C(k')=Q C(k) Q^T, not component-wise C_ii(k')=C_ii(k)",
        "reported_phase_example":diagnose_stage8n_phase_example(),
        "correct_dimensional_unit_frozen_2D_gamma":"Ang^-1 eV^-1 T^-1 before any physical charge/magnetic-source prefactors",
        "physical_gamma_computed":False,"new_vasp_pbs_jobs":0
    }

if __name__=="__main__":
    import json
    print(json.dumps(audit_summary(),indent=2))
