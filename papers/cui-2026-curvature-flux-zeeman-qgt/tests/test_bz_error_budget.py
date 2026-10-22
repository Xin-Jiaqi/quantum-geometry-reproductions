"""Round21: conditional cancellation error and multi-grid aliasing controls."""
import unittest
import numpy as np
from flux_lab.hexagonal_bz import INV_BASIS
from flux_lab.primitive_bz import primitive_cell_quadrature
from flux_lab.bz_error_budget import (
    cancellation_condition,summarize_independent_estimates,
    shifted_midpoint_ensemble,_unit_allowed)

class CancellationCondition(unittest.TestCase):
    def test_exact_condition_identity_from_triangle_inequality(self):
        # Tiny nonzero net of a much larger signed integrand amplifies errors.
        w=np.array([.5,.5])
        F=np.array([1.,-.999])
        delta=np.array([1e-5,1e-5])
        I=w@F;L=w@abs(F);r=abs(I)/L
        eps=(w@abs(delta))/L
        condition=cancellation_condition(r,desired_relative_net_error=.01)
        self.assertGreater(condition["condition_amplification"],1000)
        self.assertLessEqual(abs(w@delta)/abs(I),eps/r+1e-12)
        self.assertAlmostEqual(condition["required_relative_l1_integrand_accuracy"],.01*r)

    def test_model_specific_1percent_accuracy_requirement(self):
        six=cancellation_condition(.014248470102075281)
        mirror=cancellation_condition(.7887923373042899)
        self.assertAlmostEqual(six["required_relative_l1_integrand_accuracy"],
                               .0001424847010207528)
        self.assertAlmostEqual(mirror["required_relative_l1_integrand_accuracy"],
                               .007887923373042899)
        self.assertGreater(six["condition_amplification"],70)
        self.assertFalse(six["certified"])

    def test_exact_zero_cannot_have_relative_net_error(self):
        d=cancellation_condition(0.)
        self.assertIsNone(d["condition_amplification"])
        self.assertEqual(d["required_relative_l1_integrand_accuracy"],0.)

    def test_invalid_conditions_rejected(self):
        for r in (-.1,1.1,float("nan")):
            with self.assertRaises(ValueError):cancellation_condition(r)
        with self.assertRaises(ValueError):
            cancellation_condition(.4,desired_relative_net_error=0.)

    def test_independent_ensemble_spreads_are_complex_not_magnitude_only(self):
        # Same |response| can have opposite phase and large tensor error.
        x=summarize_independent_estimates(
            coarse_midpoint=[1+0j,-1+0j],
            fine_midpoint=[1j,-1j],
            hex_reference=1+0j,primitive_gauss=-1+0j)
        self.assertGreaterEqual(x["independent_method_spread"],2.)
        self.assertGreaterEqual(x["fine_shift_spread"],2.)
        self.assertEqual(x["numerical_convergence_status"],"UNVERIFIED")
        self.assertFalse(x["is_rigorous_error_upper_bound"])

    def test_coincident_estimates_are_not_convergence_proof(self):
        d=summarize_independent_estimates(
            coarse_midpoint=[0.,0.],fine_midpoint=[0.,0.],
            hex_reference=0.,primitive_gauss=0.)
        self.assertEqual(d["empirical_disagreement"],0.)
        self.assertEqual(d["numerical_convergence_status"],"UNVERIFIED")

    def test_periodic_aliasing_negative_control(self):
        # cos(8*q1) has exact primitive torus average 0, but an 8x8
        # unshifted midpoint grid samples a single phase (aliasing).
        k0,w0,_=primitive_cell_quadrature(8,method="midpoint",shift=(0.,0.))
        k1,w1,_=primitive_cell_quadrature(8,method="midpoint",shift=(.25,0.))
        q0=k0@INV_BASIS.T;q1=k1@INV_BASIS.T
        f0=float(np.sum(w0*np.cos(8*q0[:,0])))
        f1=float(np.sum(w1*np.cos(8*q1[:,0])))
        self.assertGreater(abs(f0),.999)
        self.assertLess(abs(f1),1e-13)
        self.assertGreater(abs(f0-f1),.999)

    def test_symmetry_allowed_and_forbidden_channels_separated(self):
        self.assertTrue(_unit_allowed("6","even",(1,1)))
        self.assertTrue(_unit_allowed("6m'm'","even",(0,1)))
        self.assertFalse(_unit_allowed("6","even",(0,2)))

    def test_small_real_Kubo_sampling_comparison_preserves_full_tensor(self):
        d=shifted_midpoint_ensemble("6",component=(1,1),parity="even",
              coarse_order=2,fine_order=3,hex_order=2,gauss_order=3,
              shifts=((0.,0.),(.25,-.2)))
        self.assertEqual(len(d["fine_midpoint"]),2)
        self.assertTrue(d["symmetry_allowed"])
        self.assertTrue(np.isfinite(d["comparison"]["empirical_disagreement"]))
        self.assertEqual(d["comparison"]["numerical_convergence_status"],"UNVERIFIED")

    def test_small_magnetic_group_sampling_does_not_project_results(self):
        d=shifted_midpoint_ensemble("6m'm'",component=(0,1),parity="even",
              coarse_order=2,fine_order=3,hex_order=2,gauss_order=3,
              shifts=((0.,0.),(.2,-.25)))
        self.assertTrue(d["symmetry_allowed"])
        self.assertGreaterEqual(d["gauss_projector_residual"],0.)
        self.assertFalse(d["comparison"]["is_rigorous_error_upper_bound"])

    def test_invalid_sampling_configuration_rejected(self):
        with self.assertRaises(ValueError):
            shifted_midpoint_ensemble("6",component=(1,1),parity="wrong")
        with self.assertRaises(ValueError):
            shifted_midpoint_ensemble("6",component=(1,1),
                       coarse_order=4,fine_order=3)
        with self.assertRaises(ValueError):
            shifted_midpoint_ensemble("6",component=(1,1),shifts=((0.,0.),))
        with self.assertRaises(ValueError):
            shifted_midpoint_ensemble("6",component=(1,1),shifts=((0.,0.),(.8,0.)))
        with self.assertRaises(ValueError):
            summarize_independent_estimates(coarse_midpoint=[1],
                fine_midpoint=[1,2],hex_reference=1,primitive_gauss=1)

if __name__=="__main__":
    unittest.main()
