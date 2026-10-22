"""v72 finite-q Bloch spin Zeeman source and charge continuity science gates."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.bulk_finite_q_spin_v72 import (
    Spec,spectral,liouville_independent,analytic_zeeman_mixed,field_difference,
    incorrect_identity_vertex,incorrect_local_velocity_vertex,_entries,_field_and_vertices,
    _dag,SX)

class BulkSpinFiniteQ(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec()
        cls.base=spectral(cls.s)
        cls.ana=analytic_zeeman_mixed(cls.s)
        cls.fd=field_difference(cls.s,step=1e-4)
        cls.liou=liouville_independent(cls.s)
        cls.g=replace(cls.s,orbital_gauge=.7)
        cls.gb=spectral(cls.g)
        cls.ga=analytic_zeeman_mixed(cls.g)

    def test_01_h_is_periodic_and_hermitian(self):
        h,_,_,_=_field_and_vertices(self.s)
        self.assertLess(np.linalg.norm(h-_dag(h)),1e-12)

    def test_02_spin_is_explicit_zeeman_not_peierls(self):
        a=_field_and_vertices(self.s)[0]
        b=_field_and_vertices(replace(self.s,field=self.s.field+.1))[0]
        self.assertLess(np.linalg.norm(b-a-.1*SX),1e-12)

    def test_03_fixed_k_electron_current_independent_B(self):
        a=_entries(self.s)[5]
        b=_entries(replace(self.s,field=self.s.field+.05))[5]
        self.assertLess(np.linalg.norm(a-b),1e-12)

    def test_04_finite_q_longitudinal_vertex_independent_B(self):
        a=_entries(self.s)[6]
        b=_entries(replace(self.s,field=self.s.field+.06))[6]
        self.assertLess(np.linalg.norm(a-b),1e-12)

    def test_05_no_q_zero_undefined_scalar_E_division(self):
        with self.assertRaises(ValueError):Spec(q_steps=0)

    def test_06_periodic_q_commensurate(self):
        self.assertAlmostEqual(self.base['q'],2*np.pi/24)

    def test_07_explicit_density_nonzero_at_finite_q(self):
        self.assertGreater(abs(self.base['density']),.003)

    def test_08_explicit_spin_mixed_transverse_nonzero(self):
        self.assertGreater(abs(self.ana['dB_E_transverse']),.003)

    def test_09_exact_global_density_fsum(self):
        self.assertLess(abs(self.base['global_density_fsum']),1e-13)

    def test_10_charge_ward_at_finite_q(self):
        self.assertLess(abs(self.base['ward_residual']),1e-12)

    def test_11_mixed_charge_ward_at_finite_q(self):
        self.assertLess(abs(self.ana['dB_ward_residual']),1e-11)

    def test_12_density_independent_dense_sylvester(self):
        self.assertLess(abs(self.base['density']-self.liou['density']),1e-12)

    def test_13_longitudinal_independent_dense_sylvester(self):
        self.assertLess(abs(self.base['longitudinal']-self.liou['longitudinal']),1e-12)

    def test_14_transverse_independent_dense_sylvester(self):
        self.assertLess(abs(self.base['transverse']-self.liou['transverse']),1e-12)

    def test_15_analytic_B_density_against_five_point(self):
        self.assertLess(abs(self.ana['dB_density']-self.fd['density']),2e-9)

    def test_16_analytic_B_longitudinal_against_five_point(self):
        self.assertLess(abs(self.ana['dB_longitudinal']-self.fd['longitudinal']),2e-9)

    def test_17_analytic_B_transverse_against_five_point(self):
        self.assertLess(abs(self.ana['dB_transverse']-self.fd['transverse']),2e-9)

    def test_18_analytic_B_E_normalization_against_five_point(self):
        self.assertLess(abs(self.ana['dB_E_transverse']-self.fd['electric_normalized_transverse']),2e-9)

    def test_19_full_k_dependent_gauge_density(self):
        self.assertLess(abs(self.gb['density']-self.base['density']),1e-12)

    def test_20_full_k_dependent_gauge_longitudinal(self):
        self.assertLess(abs(self.gb['longitudinal']-self.base['longitudinal']),1e-12)

    def test_21_full_k_dependent_gauge_transverse(self):
        self.assertLess(abs(self.gb['transverse']-self.base['transverse']),1e-12)

    def test_22_full_k_dependent_gauge_mixed_transverse(self):
        self.assertLess(abs(self.ga['dB_transverse']-self.ana['dB_transverse']),1e-12)

    def test_23_full_k_dependent_gauge_mixed_charge_ward(self):
        self.assertLess(abs(self.ga['dB_ward_residual']),1e-12)

    def test_24_fake_identity_density_vertex_breaks_covariance(self):
        wrong=incorrect_identity_vertex(self.g)
        self.assertGreater(abs(wrong['density']-self.base['density']),4e-4)

    def test_25_fake_identity_density_vertex_breaks_Ward(self):
        wrong=incorrect_identity_vertex(self.g)
        self.assertGreater(abs(wrong['ward_defect']),5e-4)

    def test_26_midpoint_velocity_fails_exact_finite_q_Ward(self):
        self.assertGreater(abs(incorrect_local_velocity_vertex(self.s)['ward_defect']),1e-5)

    def test_27_midpoint_vertex_larger_error_at_larger_q(self):
        d1=abs(incorrect_local_velocity_vertex(self.s)['ward_defect'])
        d2=abs(incorrect_local_velocity_vertex(replace(self.s,q_steps=2))['ward_defect'])
        self.assertGreater(d2,5*d1)

    def test_28_commensurate_second_momentum_Ward(self):
        s=replace(self.s,q_steps=2)
        self.assertLess(abs(spectral(s)['ward_residual']),1e-12)
        self.assertLess(abs(analytic_zeeman_mixed(s)['dB_ward_residual']),1e-11)

    def test_29_higher_B_finite_difference_step(self):
        fd=field_difference(self.s,step=2e-4)
        self.assertLess(abs(fd['transverse']-self.ana['dB_transverse']),1e-8)

    def test_30_momentum_dependent_spin_frame_Z_changes_with_k(self):
        a=_entries(self.g)
        self.assertGreater(np.linalg.norm(a[2]-a[3]),1e-2)

    def test_31_finite_q_spectral_value_frozen(self):
        self.assertLess(abs(self.base['electric_normalized_transverse']-(.006895888683487931+.0031287537300871094j)),1e-12)

    def test_32_spin_mixed_value_frozen(self):
        self.assertLess(abs(self.ana['dB_E_transverse']-(-.005924308951739088-.007974278348437747j)),1e-12)

    def test_33_invalid_eta_and_T(self):
        with self.assertRaises(ValueError):Spec(eta=0)
        with self.assertRaises(ValueError):Spec(temperature=0)

    def test_34_wrong_q_steps_rejected(self):
        with self.assertRaises(ValueError):Spec(n=24,q_steps=12)
        with self.assertRaises(ValueError):Spec(n=24,q_steps=.4)

    def test_35_invalid_magnetic_finite_difference(self):
        with self.assertRaises(ValueError):field_difference(self.s,step=0)
        with self.assertRaises(ValueError):field_difference(self.s,points=4)

    def test_36_negative_gauge_requires_actual_k_dependent_frame(self):
        with self.assertRaises(ValueError):incorrect_identity_vertex(self.s)

    def test_37_no_source_DC_claim(self):
        self.assertIn('NOT_SOURCE_PHE_DC',self.base['status'])
        self.assertIn('NO_ORBITAL_FIELD',self.ana['status'])

if __name__=='__main__':unittest.main()
