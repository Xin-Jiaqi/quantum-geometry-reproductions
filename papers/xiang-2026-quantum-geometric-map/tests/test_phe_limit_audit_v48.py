"""v48: scientific gates for original Eq11 (NOT source-SM claim)."""
import unittest
import numpy as np
from qgmap.dirac_mixed_rta_v46 import Spec
from qgmap.phe_limit_audit_v48 import eta_limit,chemical_potential_temperature,cutoff_sweep,_simple
from qgmap.dirac_source_tensor_v47 import derivative_source_check,tensor


class OriginalPHELimitV48(unittest.TestCase):
    def test_eta2_extrapolate_to_unity_and_truncation_control(self):
        r=eta_limit()
        self.assertEqual(r['scientific_status'],'PARTIAL')
        self.assertEqual(len(r['values']),5)
        self.assertLess(r['fit_eta_squared'][0]['absolute_intercept_minus_unity'],2e-6)
        self.assertLess(r['fit_eta_squared'][-1]['absolute_intercept_minus_unity'],8e-6)
        self.assertGreater(abs(r['values'][0]['ratio_to_main_eq11']-1),.003)

    def test_eta_squared_scaling_is_not_linear(self):
        r=eta_limit()
        pts=r['values']
        for p,q in zip(pts[1:-1],pts[2:]):
            a=abs(p['ratio_to_main_eq11']-1)
            b=abs(q['ratio_to_main_eq11']-1)
            self.assertTrue(3.6<a/b<4.3)

    def test_antisymmetric_eta_times_tensor_vanishes_quadratically(self):
        r=eta_limit()
        self.assertLess(max(r['antisymmetric_eta_scaled_by_eta2'])-min(r['antisymmetric_eta_scaled_by_eta2']),.001)
        self.assertLess(abs(r['values'][-1]['eta_antisymmetric']),2e-8)

    def test_mu_reversal_odd_and_tanh_at_two_temperatures(self):
        r=chemical_potential_temperature()
        self.assertEqual(len(r['results']),4)
        for x in r['results']:
            self.assertLess(abs(x['odd_parity_residual']),1e-10)
            self.assertLess(abs(x['positive_ratio']-1),.001)
            self.assertLess(abs(x['negative_ratio']-1),.001)

    def test_thermal_smoothing_is_not_absent(self):
        r=chemical_potential_temperature()
        plus=[x['positive_mu']['eta_symmetric'] for x in r['results'] if x['mu_abs']==.15]
        self.assertGreater(abs(plus[0]),abs(plus[1]))

    def test_UV_spread_bounded_at_fixed_beta_and_regulator(self):
        r=cutoff_sweep()
        self.assertEqual(len(r['values']),4)
        self.assertLess(r['max_ratio_spread'],5e-7)

    def test_source_full_Fermi_matrix_fd_still_independent(self):
        d=derivative_source_check(Spec(nr=18,nphi=32),step=2e-5)
        self.assertLess(d['max_source_y_fd_error'],2e-7)

    def test_named_tensor_not_xy_only(self):
        v=tensor(Spec(eta=.002,temperature=.03,nr=140,nphi=160),step=2e-5)
        self.assertAlmostEqual(v['symmetric'],(v['xy']['total']+v['yx']['total'])/2,delta=1e-10)
        self.assertLess(abs(v['ratio_to_main_eq11']-1),1e-4)

    def test_no_source_supplemental_claim_in_protocol(self):
        r=eta_limit()
        self.assertTrue(any('SM collision operator' in x for x in r['limitations']))
        self.assertEqual(r['scientific_status'],'PARTIAL')

    def test_fail_closed_parameter_controls(self):
        with self.assertRaises(ValueError):eta_limit(etas=(.001,.002,.003,.004))
        with self.assertRaises(ValueError):eta_limit(etas=(.01,.01,.005,.001))
        with self.assertRaises(ValueError):eta_limit(etas=(.01,.005,.001))
        with self.assertRaises(ValueError):chemical_potential_temperature(mus=(-.3,))
        with self.assertRaises(ValueError):cutoff_sweep(cutoffs=(-1.,2.))

if __name__=='__main__':unittest.main()
