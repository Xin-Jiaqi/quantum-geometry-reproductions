"""v52 independent original Fig2(b) 43microV and absolute SI unit chain physics gates."""
import math
import unittest
from dataclasses import replace
from scipy import constants
from qgmap.fig2_si_v52 import (
 Fig2Config,eq11_f,independent_eq11_raw,absolute_fig2_audit,
 finite_temperature_shape_audit,full_report,SOURCE_MAIN_SHA256,SOURCE_SM_SHA256)

class PublishedFig2SI(unittest.TestCase):
    def test_stable_eq11_matches_original_unsimplified_beta_sweep(self):
        for beta in (.025,.08,.125,.25,.5,.8):
            self.assertAlmostEqual(eq11_f(beta),independent_eq11_raw(beta),delta=2e-12)
            self.assertLess(eq11_f(beta),0)

    def test_e_squared_over_hbar_vs_e_squared_over_h_is_exactly_2pi(self):
        a=absolute_fig2_audit()
        self.assertAlmostEqual(a['physics_constants']['ratio_e2hbar_over_e2h'],2*math.pi,delta=3e-15)
        self.assertAlmostEqual(a['quantitative_comparisons']['literal_divided_by_e2h_variant'],2*math.pi,delta=3e-15)

    def test_independent_joule_and_ev_zeeman_energy_conversion(self):
        a=absolute_fig2_audit(); d=a['params']
        m_j=constants.physical_constants['Bohr magneton'][0]
        m_ev=constants.physical_constants['Bohr magneton in eV/T'][0]
        evhbar=constants.physical_constants['Planck constant over 2 pi in eV s'][0]
        pref_ev=d['tau_ps']*1e-12*d['g_factor']*m_ev*d['bx_tesla']/evhbar
        # CODATA tabulates mu_B in eV/T to fewer significant figures than in J/T.
        self.assertAlmostEqual(pref_ev,a['tau_zeeman_over_hbar'],delta=1e-8)
        self.assertAlmostEqual(a['zeeman_g_muB_B_ev'],
                               d['g_factor']*m_j*d['bx_tesla']/constants.e,delta=1e-12)

    def test_literal_eq11_fig2_plateau_is_about_11_not_1_75(self):
        a=absolute_fig2_audit(); y=a['literal_Eq11']['figure_axis_minus_1e4']
        self.assertAlmostEqual(y,11.02083611497,delta=5e-8)
        self.assertGreater(abs(y-1.75),9)

    def test_conductance_e2h_variant_predicts_published_plateau(self):
        a=absolute_fig2_audit()
        self.assertAlmostEqual(a['independent_e2h_variant']['figure_axis_minus_1e4'],
                               1.75402054470,delta=5e-8)
        self.assertLess(abs(a['quantitative_comparisons']['h_variant_plateau_minus_observed']),
                        a['published']['figure_digitization_tolerance'])

    def test_literal_voltage_independently_closes_unit_chain_and_exceeds_43(self):
        a=absolute_fig2_audit()
        literal=abs(a['literal_Eq11']['sigma_S'])
        independent_v=literal*100e-6*1e4*1e3*1e6
        self.assertAlmostEqual(a['literal_Eq11']['predicted_voltage_microV'],
                               independent_v,delta=1e-8)
        self.assertAlmostEqual(independent_v,268.26200776,delta=1e-5)
        self.assertGreater(independent_v/43,6)

    def test_e2h_variant_predicts_text_43_microvolts_without_fit(self):
        a=absolute_fig2_audit()
        self.assertAlmostEqual(a['independent_e2h_variant']['predicted_voltage_microV'],42.69522458,delta=1e-5)
        self.assertLess(abs(a['independent_e2h_variant']['predicted_voltage_microV']-43),.4)
        self.assertAlmostEqual(a['literal_Eq11']['predicted_voltage_microV']/
            a['independent_e2h_variant']['predicted_voltage_microV'],2*math.pi,delta=1e-13)

    def test_factor_2pi_persists_with_beta_and_tau_g_field_variation(self):
        for beta in (.05,.125,.3,.5):
            for tau in (.2,.5,.9):
                for g in (2,10):
                    a=absolute_fig2_audit(Fig2Config(beta=beta,tau_ps=tau,g_factor=g))
                    self.assertAlmostEqual(a['literal_Eq11']['figure_axis_minus_1e4']/
                       a['independent_e2h_variant']['figure_axis_minus_1e4'],
                       2*math.pi,delta=1e-13)

    def test_independent_parameter_variation_no_hidden_geometry_fit(self):
        a=absolute_fig2_audit()
        other=absolute_fig2_audit(replace(Fig2Config(),ey_v_per_m=2e4,lateral_um=50))
        self.assertAlmostEqual(a['literal_Eq11']['predicted_voltage_microV'],
                               other['literal_Eq11']['predicted_voltage_microV'],delta=1e-9)
        self.assertAlmostEqual(a['literal_Eq11']['figure_axis_minus_1e4'],
                               other['literal_Eq11']['figure_axis_minus_1e4'],delta=1e-14)

    def test_published_green_blue_finite_T_curves_match_tanh_shape_unfitted(self):
        d=finite_temperature_shape_audit()
        self.assertEqual(len(d['input_samples']),8)
        self.assertLess(d['max_abs_residual'],.025)
        self.assertIn('model',d['limitation'])
        self.assertEqual({x['T_K'] for x in d['input_samples']},{20,40})

    def test_thermal_width_changes_with_temperature_not_constant_fit(self):
        d=finite_temperature_shape_audit()
        vals={(x['T_K'],x['mu_meV']):x['no_amplitude_fit_prediction'] for x in d['input_samples']}
        self.assertGreater(vals[(20,5)],vals[(40,5)])
        self.assertAlmostEqual(vals[(20,-5)],-vals[(20,5)],delta=1e-14)
        self.assertAlmostEqual(vals[(20,10)],-vals[(20,-10)],delta=1e-14)

    def test_no_source_erratum_or_physical_sign_inference(self):
        z=full_report()
        self.assertEqual(z['scientific_status'],'PARTIAL')
        self.assertEqual(z['absolute_SI']['scientific_status'],'PARTIAL')
        self.assertTrue(any('erratum' in s for s in z['absolute_SI']['not_proven']))
        self.assertIn('hbar',z['absolute_SI']['published']['figure_axis_unit_printed'])
        self.assertEqual(len(SOURCE_MAIN_SHA256),64)
        self.assertEqual(len(SOURCE_SM_SHA256),64)

    def test_fail_closed_invalid_units_or_beta(self):
        for args in ({'beta':0},{'beta':1},{'g_factor':0},{'tau_ps':-.1},
                     {'lateral_um':0},{'bx_tesla':0},{'ey_v_per_m':0},{'resistance_ohm':0}):
            with self.assertRaises(ValueError):Fig2Config(**args)
        with self.assertRaises(ValueError):eq11_f(0)
        with self.assertRaises(ValueError):independent_eq11_raw(1)
        with self.assertRaises(ValueError):finite_temperature_shape_audit(samples=((0,5,100),))

if __name__=='__main__':unittest.main()
