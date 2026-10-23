"""Stage 7i physics preflight: exact spatial symmetry restrictions."""
import unittest
import numpy as np
from bsf_nmee.symmetry_dft_gate import examples,spatial_transform


class StaticMixedTensorSymmetryTests(unittest.TestCase):
    def test_fixed_stack_centrosymmetry_and_horizontal_mirror_forbid_xzx(self):
        x=examples()
        self.assertTrue(x["C1"])
        self.assertTrue(x["C3"])
        self.assertFalse(x["horizontal_mirror_mz"])
        self.assertFalse(x["inversion"])
        self.assertFalse(x["C3_plus_mz"])

    def test_time_reversal_does_not_forbid_zeeman_derivative(self):
        self.assertFalse(examples()["TR_alone_forbids"])

    def test_orthogonal_transform_of_two_axial_one_polar(self):
        M=np.zeros((3,3,3));M[0,2,0]=1
        R=np.diag([1,1,-1])
        self.assertAlmostEqual(spatial_transform(M,R)[0,2,0],-1)


if __name__=='__main__':
    unittest.main()
