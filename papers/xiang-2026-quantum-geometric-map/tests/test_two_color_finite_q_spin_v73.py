"""v73 independent finite-q two-color order/regulator and static preparation physics."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.two_color_finite_q_spin_v73 import (
    Spec, two_color, two_color_dense_independent, static_magnetic_derivative,
    regulator_order_scan, _eig_data, _source_and_currents, _spectral_resolvent,
    _population_matrix, _parent)
from qgmap.bulk_finite_q_spin_v72 import analytic_zeeman_mixed

class V73TwoColorGates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(n=24)
        cls.spectral=two_color(cls.s)
        cls.dense=two_color_dense_independent(cls.s)
        cls.gauge=two_color(replace(cls.s,orbital_gauge=.71))
        cls.static_001=static_magnetic_derivative(replace(cls.s,omega_b=0,eta_b=.001))
        cls.static_0001=static_magnetic_derivative(replace(cls.s,omega_b=0,eta_b=.0001))
        cls.static_00001=static_magnetic_derivative(replace(cls.s,omega_b=0,eta_b=.00001))
        cls.static_frame=two_color(replace(cls.s,omega_b=0,eta_b=.0001,orbital_gauge=.71),True)
        cls.static_plain=two_color(replace(cls.s,omega_b=0,eta_b=.0001),True)
        cls.scan=regulator_order_scan(cls.s)

    def test_01_distinct_nonzero_orders(self):
        a=self.spectral
        self.assertGreater(abs(a['b_then_e']['normalized_E_transverse']),.005)
        self.assertGreater(abs(a['e_then_b']['normalized_E_transverse']),.02)

    def test_02_two_orderings_sum_to_total(self):
        a=self.spectral
        for key in ('density','longitudinal','transverse','normalized_E_transverse'):
            self.assertLess(abs(a['b_then_e'][key]+a['e_then_b'][key]-a['total'][key]),1e-13)

    def test_03_full_two_color_integrated_source_zero(self):
        self.assertLess(abs(self.spectral['total_source_integrated']),1e-12)

    def test_04_full_two_color_charge_Ward(self):
        self.assertLess(abs(self.spectral['ward_residual']),1e-12)

    def test_05_omega_sum_is_two_independent_poles(self):
        a=self.spectral
        self.assertLess(abs(a['frequency_sum']-a['electric_pole']-a['magnetic_pole']),1e-15)
        self.assertAlmostEqual(a['frequency_sum'].imag,.075,places=14)

    def test_06_different_etas_are_distinct(self):
        a=two_color(replace(self.s,eta_b=.02))['total']['normalized_E_transverse']
        b=two_color(replace(self.s,eta_e=.02))['total']['normalized_E_transverse']
        self.assertGreater(abs(a-b),.002)

    def test_07_incorrect_total_pole_changes_result(self):
        a=two_color(replace(self.s,eta_b=.03))['total']['normalized_E_transverse']
        b=two_color(replace(self.s,eta_b=.002))['total']['normalized_E_transverse']
        self.assertGreater(abs(a-b),.002)

    def test_08_e_first_order_not_redundant(self):
        a=self.spectral
        self.assertGreater(abs(a['total']['normalized_E_transverse']-a['b_then_e']['normalized_E_transverse']),.02)

    def test_09_b_first_order_not_redundant(self):
        a=self.spectral
        self.assertGreater(abs(a['total']['normalized_E_transverse']-a['e_then_b']['normalized_E_transverse']),.005)

    def test_10_spectral_dense_ebranch(self):
        a=self.spectral['e_then_b'];b=self.dense['e_then_b']
        for k in a:self.assertLess(abs(a[k]-b[k]),1e-12,k)

    def test_11_spectral_dense_bebranch(self):
        a=self.spectral['b_then_e'];b=self.dense['b_then_e']
        for k in a:self.assertLess(abs(a[k]-b[k]),1e-12,k)

    def test_12_spectral_dense_total(self):
        a=self.spectral['total'];b=self.dense['total']
        for k in a:self.assertLess(abs(a[k]-b[k]),1e-12,k)

    def test_13_gauge_e_then_b(self):
        for k in self.spectral['e_then_b']:
            self.assertLess(abs(self.spectral['e_then_b'][k]-self.gauge['e_then_b'][k]),1e-12,k)

    def test_14_gauge_b_then_e(self):
        for k in self.spectral['b_then_e']:
            self.assertLess(abs(self.spectral['b_then_e'][k]-self.gauge['b_then_e'][k]),1e-12,k)

    def test_15_gauge_total(self):
        for k in self.spectral['total']:
            self.assertLess(abs(self.spectral['total'][k]-self.gauge['total'][k]),1e-12,k)

    def test_16_gauge_ward(self):
        self.assertLess(abs(self.gauge['ward_residual']),1e-12)

    def test_17_q2_mixed_ward(self):
        for qstep in (2,3):
            self.assertLess(abs(two_color(replace(self.s,q_steps=qstep))['ward_residual']),1e-12)

    def test_18_omega_b_variation_necessary(self):
        a=self.spectral['total']['normalized_E_transverse']
        b=two_color(replace(self.s,omega_b=0.))['total']['normalized_E_transverse']
        self.assertGreater(abs(a-b),.01)

    def test_19_frozen_two_color_value(self):
        z=self.spectral['total']['normalized_E_transverse']
        self.assertLess(abs(z-(-.04200892243622867-.009686339341006673j)),1e-12)

    def test_20_frozen_e_then_b_value(self):
        z=self.spectral['e_then_b']['normalized_E_transverse']
        self.assertLess(abs(z-(-.043217183392484496-.019011725385956996j)),1e-12)

    def test_21_frozen_b_then_e_value(self):
        z=self.spectral['b_then_e']['normalized_E_transverse']
        self.assertLess(abs(z-(.0012082609562558142+.009325386044950316j)),1e-12)

    def test_22_static_population_is_missing_in_collisionless(self):
        self.assertGreater(abs(self.static_00001['missing_population_signal']),.02)

    def test_23_static_collisionless_error_nonzero(self):
        self.assertGreater(self.static_00001['collisionless_error'],.02)

    def test_24_prepared_static_improves_vs_collisionless(self):
        for a in (self.static_001,self.static_0001,self.static_00001):
            self.assertLess(a['prepared_error'],.1*a['collisionless_error'])

    def test_25_prepared_static_error_scales_eta_b(self):
        self.assertLess(self.static_0001['prepared_error'],self.static_001['prepared_error']*.13)
        self.assertLess(self.static_00001['prepared_error'],self.static_0001['prepared_error']*.13)

    def test_26_prepared_static_reaches_analytic_B_derivative(self):
        self.assertLess(self.static_00001['prepared_error'],3e-6)

    def test_27_static_prepared_target_v72(self):
        a=self.static_00001['exact_prepared_dB']
        b=analytic_zeeman_mixed(_parent(self.s))['dB_E_transverse']
        self.assertLess(abs(a-b),1e-12)

    def test_28_prepared_static_field_frame(self):
        for key in ('b_then_e','e_then_b','total'):
            for name in ('density','longitudinal','transverse','normalized_E_transverse'):
                self.assertLess(abs(self.static_frame[key][name]-self.static_plain[key][name]),1e-12)

    def test_29_prepared_static_ward(self):
        self.assertLess(abs(self.static_plain['ward_residual']),1e-12)

    def test_30_dense_prepared_static(self):
        s=replace(self.s,omega_b=0.,eta_b=.0001)
        spectral=two_color(s,True)
        dense=two_color_dense_independent(s,True)
        for key in ('b_then_e','e_then_b','total'):
            self.assertLess(abs(spectral[key]['normalized_E_transverse']-
                                dense[key]['normalized_E_transverse']),1e-12)

    def test_31_both_eta_poles_for_static_preparation(self):
        s=replace(self.s,omega_b=0.,eta_b=.001)
        a=two_color(s,True)
        self.assertAlmostEqual(a['frequency_sum'].imag,.046,places=14)

    def test_32_no_fake_static_equilibrium_at_nonzero_b_frequency(self):
        with self.assertRaises(ValueError):two_color(self.s,True)
        with self.assertRaises(ValueError):two_color_dense_independent(self.s,True)
        with self.assertRaises(ValueError):static_magnetic_derivative(self.s)

    def test_33_static_prepare_flag_exposed(self):
        self.assertFalse(self.spectral['prepared_static_population'])
        self.assertTrue(self.static_plain['prepared_static_population'])

    def test_34_no_source_tau_pole_static_claim(self):
        self.assertIn('NOT_SOURCE_TAU_POLE',self.static_plain['status'])
        self.assertIn('NOT_SOURCE_PHE_DC',self.spectral['status'])

    def test_35_regulator_scan_does_not_become_source_tau_at_fixed_q(self):
        self.assertIn('NO_SOURCE_TAU_IDENTIFICATION',self.scan['status'])
        z=[abs(v['eta_e_times_value']) for v in self.scan['values'].values()]
        self.assertLess(z[-1],.5*z[0])

    def test_36_regulator_scan_values_finite(self):
        for d in self.scan['values'].values():
            self.assertTrue(np.isfinite([d['value'].real,d['value'].imag]).all())

    def test_37_finite_q_E_normalization(self):
        self.assertAlmostEqual(self.spectral['q'],2*np.pi/24)
        self.assertLess(abs(self.spectral['total']['transverse']/(-1j*self.spectral['q'])-
                            self.spectral['total']['normalized_E_transverse']),1e-15)

    def test_38_finite_q_q0_rejected(self):
        with self.assertRaises(ValueError):Spec(q_steps=0)

    def test_39_etas_positive(self):
        for kw in ({'eta_e':0.},{'eta_b':0.},{'temperature':-1.}):
            with self.assertRaises(ValueError):Spec(**kw)

    def test_40_invalid_q_steps(self):
        for qstep in (-1,12,1.1):
            with self.assertRaises(ValueError):Spec(q_steps=qstep)

    def test_41_finite_all_params(self):
        with self.assertRaises(ValueError):Spec(omega_b=np.inf)
        with self.assertRaises(ValueError):Spec(eta_e=np.nan)

    def test_42_distinct_omega_e_and_omega_b(self):
        a=self.spectral
        self.assertAlmostEqual(a['electric_pole'].real,.27)
        self.assertAlmostEqual(a['magnetic_pole'].real,.12)
        self.assertAlmostEqual(a['frequency_sum'].real,.39)

    def test_43_independent_current_vertex_at_field_zero(self):
        a=two_color(replace(self.s,B0=0.))['total']['normalized_E_transverse']
        b=self.spectral['total']['normalized_E_transverse']
        self.assertGreater(abs(a-b),.0008)

    def test_44_output_regulator_different_from_input(self):
        s=replace(self.s,eta_e=.038,eta_b=.016)
        a=two_color(s)
        self.assertAlmostEqual(a['frequency_sum'].imag,.054)

    def test_45_two_color_total_does_not_equal_static_v72_derivative(self):
        a=self.spectral['total']['normalized_E_transverse']
        b=analytic_zeeman_mixed(_parent(self.s))['dB_E_transverse']
        self.assertGreater(abs(a-b),.01)

if __name__=='__main__':unittest.main()
