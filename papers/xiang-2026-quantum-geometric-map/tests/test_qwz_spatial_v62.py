"""v62 independent science acceptance: spatial NEGF and analytic-vs-Kubo bulk charge."""
import unittest
import numpy as np
from qgmap.qwz_spatial_v62 import qwz_width_leads,spatial_charge_response
from qgmap.qwz_device_v61 import qwz_four_terminal,TX,TY
from qgmap.qwz_charge_kubo_v62 import qwz_charge_kubo
from qgmap.periodic_kubo_v59 import periodic_chern

class V62SourceIndependent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.top=spatial_charge_response(7,-1,3)
        cls.bottom=spatial_charge_response(7,1,3)
        cls.trivial=spatial_charge_response(5,3,3)
        cls.full=spatial_charge_response(5,-1,5)

    def test_01_same_full_width_v61(self):
        v61=qwz_four_terminal(5,-1,eta=1e-4)
        self.assertLess(np.max(abs(self.full['transmissions']-v61['transmissions'])),2e-12)
        self.assertAlmostEqual(self.full['Rxy_h_over_q2'],v61['hall_resistance_h_over_q2'],places=11)

    def test_02_real_narrow_principal_hoppings(self):
        lead=qwz_width_leads(7,7,-1,3)
        self.assertEqual(lead[0].h00.shape,(6,6))
        self.assertTrue(np.allclose(lead[0].h00[:2,2:4],TY))
        self.assertTrue(np.allclose(lead[1].h00[:2,2:4],TX))
        self.assertTrue(np.allclose(lead[0].h01[:2,:2],TX.conj().T))
        self.assertTrue(np.allclose(lead[1].h01[:2,:2],TY))

    def test_03_contiguous_contact_support(self):
        leads=qwz_width_leads(7,7,-1,3)
        left=np.nonzero(np.any(abs(leads[0].coupling)>0,axis=1))[0]
        self.assertEqual(left.tolist(),[2*14,2*14+1,2*21,2*21+1,2*28,2*28+1])

    def test_04_invalid_width_refused(self):
        for w in (1,8,1.5):
            with self.assertRaises(ValueError):qwz_width_leads(7,7,-1,w)

    def test_05_invalid_offset_refused(self):
        with self.assertRaises(ValueError):qwz_width_leads(7,7,-1,3,offsets=(0,3,0,0))
        with self.assertRaises(ValueError):qwz_width_leads(7,7,-1,3,offsets=(0,0,0))

    def test_06_invalid_eta_refused(self):
        with self.assertRaises(ValueError):qwz_width_leads(7,7,-1,3,eta=0)

    def test_07_exact_spectral_sum_identity(self):
        self.assertLess(self.top['spectral_identity_error'],1e-10)

    def test_08_positive_ldos(self):
        self.assertGreaterEqual(self.top['spectral_ldos'].min(),-1e-12)

    def test_09_edge_larger_than_center_ldos(self):
        d=self.top['spectral_ldos']
        self.assertGreater(d[0,0],100*d[3,3])

    def test_10_interior_continuity(self):
        self.assertLess(np.max(abs(self.top['bond_divergence'][1:-1,1:-1])),1e-10)

    def test_11_total_bond_divergence(self):
        self.assertLess(abs(np.sum(self.top['bond_divergence'])),1e-10)

    def test_12_probe_currents_zero(self):
        self.assertLess(self.top['floating_current_error'],1e-12)
        self.assertLess(self.top['total_lead_current_error'],1e-12)

    def test_13_far_from_side_probes_cut_equals_drive(self):
        self.assertLess(abs(self.top['x_cut_currents'][0]-self.top['drive_current_q2_over_h']),1e-10)

    def test_14_invasive_probe_current_redistribution(self):
        d=self.top
        self.assertGreater(abs(d['x_cut_currents'][0]-d['x_cut_currents'][3]),.8)
        self.assertLess(abs(d['x_cut_currents'][3]),.05)
        self.assertLess(abs(d['Rxy_h_over_q2']-1),.003)

    def test_15_opposite_chiral_sign(self):
        self.assertLess(abs(self.top['Rxy_h_over_q2']+self.bottom['Rxy_h_over_q2']),1e-10)
        self.assertLess(abs(self.top['drive_current_q2_over_h']-self.bottom['drive_current_q2_over_h']),1e-10)

    def test_16_trivial_current_guard(self):
        self.assertIsNone(self.trivial['Rxy_h_over_q2'])
        self.assertLess(abs(self.trivial['drive_current_q2_over_h']),1e-7)

    def test_17_contact_position_changes_local_flow(self):
        changed=spatial_charge_response(7,-1,3,offsets=(0,1,0,-1))
        self.assertGreater(abs(changed['x_cut_currents'][3]-self.top['x_cut_currents'][3]),.3)
        self.assertLess(abs(changed['Rxy_h_over_q2']-1),.003)

    def test_18_local_gauge_and_contact_covariance(self):
        d=spatial_charge_response(5,-1,3)
        phase=np.linspace(-1.,1.,2*5*5)
        g=spatial_charge_response(5,-1,3,orbital_phases=phase)
        self.assertLess(np.max(abs(d['Jx']-g['Jx'])),1e-10)
        self.assertLess(np.max(abs(d['Jy']-g['Jy'])),1e-10)
        self.assertLess(np.max(abs(d['spectral_ldos']-g['spectral_ldos'])),1e-10)
        self.assertAlmostEqual(d['Rxy_h_over_q2'],g['Rxy_h_over_q2'],places=11)

    def test_19_malformed_local_gauge(self):
        with self.assertRaises(ValueError):spatial_charge_response(5,-1,3,orbital_phases=np.arange(4))

    def test_20_negative_metallic_plateau_claim(self):
        d=spatial_charge_response(5,-1,3,mu=1.5)
        c=qwz_charge_kubo(-1,1.5,0,31)['sigma_interband_q2_over_h']
        self.assertLess(abs(d['Rxy_h_over_q2']),.2)
        self.assertGreater(abs(c),.04)
        self.assertGreater(abs(d['Rxy_h_over_q2']-1/c),1)

    def test_21_qwz_two_independent_kubo_oracles(self):
        for m in (-1.,1.,3.):
            d=qwz_charge_kubo(m,.2,.16,21)
            self.assertLess(abs(d['sigma_pauli_q2_over_h']-d['sigma_interband_q2_over_h']),1e-12)

    def test_22_signed_periodic_quantization(self):
        for m,sign in ((-1,1),(1,-1),(3,0)):
            k=qwz_charge_kubo(m,0,0,31)
            self.assertLess(abs(k['sigma_pauli_q2_over_h']-sign),1e-5)
            self.assertLess(abs(k['sigma_interband_q2_over_h']-periodic_chern(m,31)['chern_spectral']),1e-11)

    def test_23_temperature_suppression(self):
        c=qwz_charge_kubo(-1,0,0,31)['sigma_pauli_q2_over_h']
        warm=qwz_charge_kubo(-1,0,1.,31)['sigma_pauli_q2_over_h']
        self.assertGreater(c,.999)
        self.assertLess(warm,.6)

    def test_24_mu_particle_hole_symmetry(self):
        a=qwz_charge_kubo(-1,-1.5,0,31)
        b=qwz_charge_kubo(-1,1.5,0,31)
        self.assertLess(abs(a['sigma_pauli_q2_over_h']-b['sigma_pauli_q2_over_h']),1e-12)

    def test_25_bulk_metal_noninteger(self):
        d=qwz_charge_kubo(-1,1.5,0,31)
        self.assertLess(abs(d['sigma_pauli_q2_over_h']),.2)
        self.assertGreater(abs(d['sigma_pauli_q2_over_h']),.01)

    def test_26_thermal_mu_combined(self):
        d=qwz_charge_kubo(-1,.35,.25,21)
        self.assertLess(abs(d['sigma_pauli_q2_over_h']-d['sigma_interband_q2_over_h']),1e-12)

    def test_27_critical_gap_refused(self):
        for m in (-2.,0.,2.):
            with self.assertRaises(ValueError):qwz_charge_kubo(m)

    def test_28_invalid_temperature_refused(self):
        with self.assertRaises(ValueError):qwz_charge_kubo(-1,0,-.1)

    def test_29_invalid_mesh_refused(self):
        with self.assertRaises(ValueError):qwz_charge_kubo(-1,0,0,8)

    def test_30_eta_robust_narrow_device(self):
        a=spatial_charge_response(5,-1,3,eta=1e-3)
        b=spatial_charge_response(5,-1,3,eta=1e-4)
        self.assertLess(abs(b['Rxy_h_over_q2']-1),.001)
        self.assertLess(abs(a['Rxy_h_over_q2']-1),.001)

if __name__=='__main__':unittest.main()
