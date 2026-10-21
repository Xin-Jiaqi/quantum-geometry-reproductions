"""v63 independent acceptance for non-invasive QWZ voltage probes."""
import unittest
import numpy as np
from qgmap.qwz_probe_limit_v63 import probe_strength_response
from qgmap.qwz_spatial_v62 import spatial_charge_response
from qgmap.qwz_charge_kubo_v62 import qwz_charge_kubo


def cut_deviation(d):
    return float(np.max(np.abs(d["x_cut_currents"]-d["drive_current_q2_over_h"])))


class TestV63ProbeLimit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s1 = probe_strength_response(7,-1,3,probe_scale=1.0,probe_width=3)
        cls.s03 = probe_strength_response(7,-1,3,probe_scale=.3,probe_width=3)
        cls.s01 = probe_strength_response(7,-1,3,probe_scale=.1,probe_width=3)
        cls.s003 = probe_strength_response(7,-1,3,probe_scale=.03,probe_width=3)
        cls.s001 = probe_strength_response(7,-1,3,probe_scale=.01,probe_width=3)
        cls.s0001 = probe_strength_response(7,-1,3,probe_scale=.001,probe_width=3)
        cls.s00001 = probe_strength_response(7,-1,3,probe_scale=.0001,probe_width=3)
        cls.v62 = spatial_charge_response(7,-1,3)
        cls.full_v63 = probe_strength_response(5,-1,5,probe_scale=1.,probe_width=5)
        cls.full_v62 = spatial_charge_response(5,-1,5)
        cls.eta_hi = probe_strength_response(7,-1,3,eta=1e-3,probe_scale=.001)
        cls.eta_mid = probe_strength_response(7,-1,3,eta=1e-4,probe_scale=.001)
        cls.eta_low = probe_strength_response(7,-1,3,eta=1e-5,probe_scale=.001)
        cls.sizes = [probe_strength_response(L,-1,3,eta=1e-4,probe_scale=1e-4)
                     for L in (5,7,9,11)]
        cls.widths_finite = [probe_strength_response(7,-1,3,probe_scale=1.,probe_width=w)
                             for w in (2,3,5,7)]
        cls.widths_weak = [probe_strength_response(7,-1,3,probe_scale=1e-4,probe_width=w)
                           for w in (2,3,5,7)]
        cls.shifted = probe_strength_response(7,-1,3,probe_scale=1e-4,offsets=(0,1,0,-1))
        cls.reverse = probe_strength_response(7,1,3,probe_scale=1e-4)
        cls.trivial = probe_strength_response(7,3,3,probe_scale=1e-4)
        cls.metal_device = probe_strength_response(5,-1,3,mu=1.5,probe_scale=1e-4)
        cls.metal_bulk = qwz_charge_kubo(-1,1.5,0.,31)

    def test_01_unit_scale_reproduces_v62_default(self):
        self.assertLess(np.max(abs(self.s1["transmissions"]-self.v62["transmissions"])),2e-12)
        self.assertAlmostEqual(self.s1["Rxy_h_over_q2"],self.v62["Rxy_h_over_q2"],places=12)

    def test_02_full_width_v61_limit_preserved(self):
        self.assertLess(np.max(abs(self.full_v63["transmissions"]-self.full_v62["transmissions"])),2e-12)
        self.assertAlmostEqual(self.full_v63["Rxy_h_over_q2"],self.full_v62["Rxy_h_over_q2"],places=12)

    def test_03_drive_contact_amplitudes_fixed(self):
        self.assertTrue(np.allclose(self.s1["coupling_norms"][[0,2]],
                                    self.s00001["coupling_norms"][[0,2]],rtol=0,atol=1e-14))

    def test_04_probe_amplitude_scales_linearly(self):
        ratio=self.s01["coupling_norms"][[1,3]]/self.s1["coupling_norms"][[1,3]]
        self.assertTrue(np.allclose(ratio,[.1,.1],rtol=1e-13,atol=1e-14))

    def test_05_probe_broadening_scales_quadratically(self):
        ratio=self.s01["gamma_norms"][[1,3]]/self.s1["gamma_norms"][[1,3]]
        self.assertTrue(np.allclose(ratio,[.01,.01],rtol=2e-8,atol=1e-12))

    def test_06_probe_selfenergy_scales_quadratically(self):
        ratio=self.s01["selfenergy_norms"][[1,3]]/self.s1["selfenergy_norms"][[1,3]]
        self.assertTrue(np.allclose(ratio,[.01,.01],rtol=2e-8,atol=1e-12))

    def test_07_drive_broadening_is_invariant(self):
        self.assertTrue(np.allclose(self.s1["gamma_norms"][[0,2]],
                                    self.s00001["gamma_norms"][[0,2]],rtol=0,atol=1e-12))

    def test_08_probe_scale_zero_is_not_a_measurement(self):
        with self.assertRaises(ValueError): probe_strength_response(7,-1,3,probe_scale=0.)

    def test_09_negative_probe_scale_refused(self):
        with self.assertRaises(ValueError): probe_strength_response(7,-1,3,probe_scale=-.1)

    def test_10_nonfinite_probe_scale_refused(self):
        with self.assertRaises(ValueError): probe_strength_response(7,-1,3,probe_scale=float("nan"))

    def test_11_invalid_probe_width_refused(self):
        for w in (1,8,2.5):
            with self.assertRaises(ValueError): probe_strength_response(7,-1,3,probe_width=w)

    def test_12_invalid_probe_offset_refused(self):
        with self.assertRaises(ValueError):
            probe_strength_response(7,-1,3,probe_width=3,offsets=(0,3,0,0))

    def test_13_current_profile_becomes_uniform_monotonically(self):
        dev=[cut_deviation(d) for d in (self.s1,self.s03,self.s01,self.s003,self.s001,self.s0001)]
        self.assertTrue(all(b<a for a,b in zip(dev,dev[1:])),dev)

    def test_14_weak_probe_current_is_uniform(self):
        self.assertLess(cut_deviation(self.s00001),1e-6)

    def test_15_hall_resistance_has_a_stable_weak_probe_limit(self):
        self.assertLess(abs(self.s00001["Rxy_h_over_q2"]-self.s0001["Rxy_h_over_q2"]),5e-8)

    def test_16_drive_current_has_a_stable_weak_probe_limit(self):
        self.assertLess(abs(self.s00001["drive_current_q2_over_h"]-
                            self.s0001["drive_current_q2_over_h"]),5e-8)

    def test_17_floating_probe_voltages_converge(self):
        self.assertLess(np.max(abs(self.s00001["voltages"]-self.s0001["voltages"])),1e-6)

    def test_18_floating_probe_net_currents_remain_zero(self):
        for d in (self.s1,self.s01,self.s00001):
            self.assertLess(d["floating_current_error"],1e-12)

    def test_19_total_lead_current_is_conserved(self):
        for d in (self.s1,self.s01,self.s00001):
            self.assertLess(d["total_lead_current_error"],1e-12)

    def test_20_interior_bond_current_continuity(self):
        for d in (self.s1,self.s01,self.s00001):
            self.assertLess(d["max_interior_divergence"],2e-14)

    def test_21_spectral_identity_holds(self):
        for d in (self.s1,self.s01,self.s00001):
            self.assertLess(d["spectral_identity_error"],2e-12)

    def test_22_ldos_is_nonnegative(self):
        self.assertGreaterEqual(float(self.s00001["spectral_ldos"].min()),-1e-12)

    def test_23_eta_regulator_sensitivity_is_bounded(self):
        vals=[d["Rxy_h_over_q2"] for d in (self.eta_hi,self.eta_mid,self.eta_low)]
        self.assertLess(max(vals)-min(vals),1.5e-5,vals)

    def test_24_finite_size_sequence_is_not_claimed_exactly_quantized(self):
        r=[d["Rxy_h_over_q2"] for d in self.sizes]
        diffs=np.diff(r)
        self.assertTrue(np.all(diffs>0),r)
        self.assertGreater(diffs[0],diffs[1])
        self.assertGreater(diffs[1],diffs[2])
        self.assertGreater(r[-1]-1,0.0029)

    def test_25_probe_width_is_independent_of_drive_width(self):
        self.assertTrue(all(d["coupling_norms"][0]==self.widths_finite[0]["coupling_norms"][0]
                            for d in self.widths_finite))
        self.assertTrue(all(d["coupling_norms"][2]==self.widths_finite[0]["coupling_norms"][2]
                            for d in self.widths_finite))

    def test_26_probe_width_changes_finite_coupling_perturbation(self):
        dev=[cut_deviation(d) for d in self.widths_finite]
        self.assertTrue(all(b>a for a,b in zip(dev,dev[1:])),dev)

    def test_27_probe_width_effect_vanishes_with_coupling(self):
        rs=[d["Rxy_h_over_q2"] for d in self.widths_weak]
        self.assertLess(max(rs)-min(rs),1e-8)
        self.assertTrue(all(cut_deviation(d)<1e-6 for d in self.widths_weak))

    def test_28_probe_position_remains_an_observable_control(self):
        self.assertGreater(abs(self.shifted["Rxy_h_over_q2"]-
                               self.s00001["Rxy_h_over_q2"]),1e-4)
        self.assertLess(cut_deviation(self.shifted),1e-6)

    def test_29_opposite_chirality_reverses_hall_sign(self):
        self.assertLess(abs(self.s00001["Rxy_h_over_q2"]+
                            self.reverse["Rxy_h_over_q2"]),1e-10)

    def test_30_trivial_insulator_fails_closed(self):
        self.assertIsNone(self.trivial["Rxy_h_over_q2"])
        self.assertFalse(self.trivial["Rxy_reportable"])
        self.assertLess(abs(self.trivial["drive_current_q2_over_h"]),1e-7)

    def test_31_metal_is_not_bulk_reciprocal_hall(self):
        sigma=self.metal_bulk["sigma_interband_q2_over_h"]
        self.assertGreater(abs(self.metal_device["Rxy_h_over_q2"]-1/sigma),1.)
        self.assertLess(abs(sigma),.1)

    def test_32_weak_limit_keeps_nonzero_probe_resolution(self):
        self.assertGreater(self.s00001["gamma_norms"][1],0.)
        self.assertGreater(self.s00001["gamma_norms"][3],0.)
        self.assertLess(self.s00001["gamma_norms"][1],1e-6)


if __name__ == "__main__":
    unittest.main()
