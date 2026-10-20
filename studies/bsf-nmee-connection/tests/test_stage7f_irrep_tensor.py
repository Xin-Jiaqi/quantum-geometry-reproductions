"""Stage7f symmetry and representation-resolved source response regression."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.irrep_channel_tensor import (
    analytic_tensor_ABC,rotate_slide,tensor_covariance_residual)


def spec_u(u):
    return replace(SlidingSpec(),shift=(.23*u,.41*u))


class Stage7FIrrepTensorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.u=.1
        cls.spec=spec_u(cls.u)
        cls.reference=analytic_tensor_ABC(cls.spec,nr=12,nt=36)

    def test_b_traces_positive_c_nematic_negative(self):
        t=self.reference['tensors']
        self.assertGreater(t['B']['trace'],0)
        self.assertLess(t['C']['nematic'][0][0],0)
        self.assertGreater(t['B']['scalar'],20*t['B']['nematic_norm'])
        self.assertGreater(t['C']['nematic_norm'],20*abs(t['C']['scalar']))
        self.assertGreater(t['A']['nematic_norm'],20*abs(t['A']['scalar']))

    def test_trace_cubic_nematic_linear_selection(self):
        below=analytic_tensor_ABC(spec_u(.02),nr=12,nt=36)['tensors']
        above=analytic_tensor_ABC(spec_u(.04),nr=12,nt=36)['tensors']
        for ch in ('A','B','C'):
            trace_ratio=above[ch]['trace']/below[ch]['trace']
            nem_ratio=above[ch]['nematic_norm']/below[ch]['nematic_norm']
            self.assertLess(abs(trace_ratio-8),.20)
            self.assertLess(abs(nem_ratio-2),.03)

    def test_all_channels_c3_and_mirror_covariance(self):
        rotated=analytic_tensor_ABC(
            rotate_slide(self.spec,120),nr=12,nt=36)
        flipped=analytic_tensor_ABC(
            replace(self.spec,shift=tuple(-np.array(self.spec.shift))),
            nr=12,nt=36)
        self.assertLess(tensor_covariance_residual(
            self.reference,rotated,120),1e-9)
        self.assertLess(tensor_covariance_residual(
            self.reference,flipped,0,odd=True),1e-9)

    def test_60_degree_trace_sign_and_full_covariance(self):
        rotated=analytic_tensor_ABC(
            rotate_slide(self.spec,60),nr=12,nt=36)
        for ch in ('A','B','C'):
            self.assertAlmostEqual(
                rotated['tensors'][ch]['trace'],
                -self.reference['tensors'][ch]['trace'],places=9)
        # Since 60deg = 240deg + inversion, response is odd under inversion.
        self.assertLess(tensor_covariance_residual(
            self.reference,rotated,240,odd=True),1e-9)


if __name__=='__main__':
    unittest.main()
