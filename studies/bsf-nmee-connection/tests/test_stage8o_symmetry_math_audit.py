"""Stage8o corrective tests for wrongly identified Stage8n SPN symmetry root cause."""
import importlib.util
import unittest
import numpy as np
from pathlib import Path

SRC=Path(__file__).resolve().parents[1]/"hpc-stage8i"/"stage8o_symmetry_math_audit.py"
spec=importlib.util.spec_from_file_location("stage8o_math_audit",SRC)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class Stage8oGeometryAndNumerics(unittest.TestCase):
    def test_claimed_cartesian_mirror_is_not_hexagonal_lattice_operation(self):
        r=m.is_crystal_lattice_rotation(m.OLD_CART_MIRROR)
        self.assertFalse(r["valid"])
        self.assertGreater(r["max_noninteger_residual"],.1)

    def test_true_hex_crystal_mirror_has_all_441_exact_k_images(self):
        R=m.cartesian_action_of_fractional_direct(m.VALID_FRAC_MIRROR)
        self.assertTrue(m.is_crystal_lattice_rotation(R)["valid"])
        self.assertTrue(np.allclose(R.T@R,np.eye(2),atol=1e-12))
        images=m.discrete_mesh_permutation(m.VALID_FRAC_MIRROR,21)
        self.assertEqual(len(images),441)
        self.assertEqual(len(set(images.values())),441)
        self.assertEqual(images[(1,3)],(3,1))

    def test_vertical_reflection_rotates_axial_spin_components(self):
        R=m.OLD_CART_MIRROR
        Q=m.axial_spin_rotation(R)
        self.assertTrue(np.allclose(Q@Q.T,np.eye(3)))
        Sx=np.diag([1.,2.]).astype(complex)
        Sy=np.diag([3.,4.]).astype(complex)
        Sz=np.diag([5.,6.]).astype(complex)
        C=m.spin_gram([Sx,Sy,Sz])
        Cnew=m.rotated_spin_gram(C,R)
        # A reflecting mirror swaps the x- and y-axial spin components.
        self.assertAlmostEqual(Cnew[0,0],C[1,1])
        self.assertAlmostEqual(Cnew[1,1],C[0,0])
        self.assertNotAlmostEqual(Cnew[0,0],C[0,0])
        self.assertAlmostEqual(np.trace(Cnew),np.trace(C))

    def test_horizontal_glide_does_not_remap_inplane_k(self):
        self.assertEqual(len(m.discrete_mesh_permutation(np.eye(2,dtype=int))),441)
        self.assertEqual(m.fractional_k_map(1/3,2/3,np.eye(2,dtype=int)),
                         (1/3,2/3))

    def test_stage8n_claimed_raw_phase_was_already_divided_by_dk_dl(self):
        p=m.diagnose_stage8n_phase_example()
        self.assertGreater(p["dk_Ang_inv"],.05)
        self.assertLess(p["dk_Ang_inv"],.2)
        self.assertLess(p["raw_loop_phase_difference_abs"],1e-10)
        self.assertGreater(p["raw_loop_phase_difference_abs"],1e-13)
        self.assertLess(p["raw_difference_over_machine_epsilon"],1e6)
        self.assertTrue(p["original_report_uses_divided_curvature_as_undivided_phase"])

    def test_correct_dimension_and_zero_extra_hpc_jobs(self):
        r=m.audit_summary()
        self.assertEqual(r["correct_dimensional_unit_frozen_2D_gamma"],"Ang^-1 eV^-1 T^-1 before any physical charge/magnetic-source prefactors")
        self.assertEqual(r["new_vasp_pbs_jobs"],0)
        self.assertFalse(r["physical_gamma_computed"])

if __name__=="__main__":
    unittest.main()
