"""v64 independent two-terminal oracle and falsifiable ordered-limit tests."""
import unittest
import numpy as np
from qgmap.qwz_virtual_probe_v64 import virtual_probe_limit
from qgmap.qwz_probe_limit_v63 import probe_strength_response
from qgmap.qwz_spatial_v62 import qwz_width_leads
from qgmap.matrix_leads_v60 import selfenergy_matrix
from qgmap.qwz_finite_marker_v60 import qwz_open


class V64OrderedLimits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref=virtual_probe_limit(7, -1., 3, 3)
        cls.small=probe_strength_response(7,probe_scale=1e-4)
        cls.medium=probe_strength_response(7,probe_scale=1e-3)
        cls.metal=virtual_probe_limit(5,mu=1.5)

    def test_01_virtual_ratio_is_signed_and_bounded(self):
        self.assertLess(abs(self.ref['Rxy_lambda_to_0_h_over_q2']-1.00215746035857),1e-9)

    def test_02_independent_probe_voltage_formula(self):
        w=self.ref['probe_spectral_weights']
        v=(w[:,0]-w[:,1])/(2*w.sum(axis=1))
        self.assertTrue(np.allclose(v,self.ref['floating_voltages_lambda_to_0'],atol=1e-13))
        self.assertLess(abs(v[0]-v[1]-self.ref['Rxy_lambda_to_0_h_over_q2']*self.ref['two_terminal_T_lr']),1e-11)

    def test_03_no_attached_probes_in_green(self):
        # Rebuild independently only two lead selfenergies, exactly the stated physics.
        H,_=qwz_open(7,7,-1.)
        leads=qwz_width_leads(7,7,-1.,3)
        sL=selfenergy_matrix(0.,leads[0])[0]
        sR=selfenergy_matrix(0.,leads[2])[0]
        gl=selfenergy_matrix(0.,leads[0])[2]
        gr=selfenergy_matrix(0.,leads[2])[2]
        g=np.linalg.inv(-H-sL-sR)
        t=float(np.trace(gl@g@gr@g.conj().T).real)
        self.assertLess(abs(t-self.ref['two_terminal_T_lr']),1e-11)

    def test_04_spectral_sum_exact(self):
        self.assertLess(self.ref['spectral_identity_error'],2e-13)

    def test_05_microscopic_left_right_detailed_balance(self):
        self.assertLess(abs(self.ref['two_terminal_T_lr']-self.ref['two_terminal_T_rl']),1e-10)

    def test_06_probe_weights_positive(self):
        self.assertGreater(self.ref['probe_spectral_weights'].min(),0)

    def test_07_virtual_top_and_bottom_weights_nonidentical(self):
        w=self.ref['probe_spectral_weights']
        self.assertGreater(abs(w[0,0]-w[1,0]),.01)

    def test_08_converges_finite_lambda(self):
        self.assertLess(abs(self.small['Rxy_h_over_q2']-self.ref['Rxy_lambda_to_0_h_over_q2']),2e-8)
        self.assertLess(abs(self.medium['Rxy_h_over_q2']-self.ref['Rxy_lambda_to_0_h_over_q2']),2e-7)

    def test_09_quadratic_voltage_convergence(self):
        a=probe_strength_response(7,probe_scale=.01)
        b=probe_strength_response(7,probe_scale=.001)
        ref=self.ref['Rxy_lambda_to_0_h_over_q2']
        ratio=(a['Rxy_h_over_q2']-ref)/(b['Rxy_h_over_q2']-ref)
        self.assertGreater(ratio,95)
        self.assertLess(ratio,105)

    def test_10_quadratic_drive_current_convergence(self):
        a=probe_strength_response(7,probe_scale=.01)
        b=probe_strength_response(7,probe_scale=.001)
        tref=self.ref['two_terminal_T_lr']
        ratio=(a['drive_current_q2_over_h']-tref)/(b['drive_current_q2_over_h']-tref)
        self.assertGreater(ratio,95)
        self.assertLess(ratio,105)

    def test_11_quadratic_internal_invasiveness(self):
        a=probe_strength_response(7,probe_scale=.01)
        b=probe_strength_response(7,probe_scale=.001)
        ea=np.max(abs(a['x_cut_currents']-a['drive_current_q2_over_h']))
        eb=np.max(abs(b['x_cut_currents']-b['drive_current_q2_over_h']))
        self.assertGreater(ea/eb,95)
        self.assertLess(ea/eb,105)

    def test_12_finite_lambda_differs_at_strong_coupling(self):
        d=probe_strength_response(7,probe_scale=1.)
        self.assertGreater(abs(d['Rxy_h_over_q2']-self.ref['Rxy_lambda_to_0_h_over_q2']),.002)
        self.assertGreater(np.max(abs(d['x_cut_currents']-d['drive_current_q2_over_h'])),.9)

    def test_13_weak_limit_probe_position_dependence(self):
        p=virtual_probe_limit(7,offsets=(0,1,0,-1))
        m=virtual_probe_limit(7,offsets=(0,-1,0,1))
        self.assertGreater(abs(p['Rxy_lambda_to_0_h_over_q2']-m['Rxy_lambda_to_0_h_over_q2']),4e-4)
        self.assertAlmostEqual(p['two_terminal_T_lr'],m['two_terminal_T_lr'],places=11)

    def test_14_strict_size_is_not_a_convergence_proof(self):
        x=[virtual_probe_limit(L)['Rxy_lambda_to_0_h_over_q2'] for L in (5,7,9,11)]
        self.assertGreater(max(x)-min(x),.001)
        self.assertTrue(all(np.isfinite(x)))

    def test_15_regulator_sweep_finite(self):
        x=[virtual_probe_limit(7,eta=e)['Rxy_lambda_to_0_h_over_q2'] for e in (1e-3,1e-4,1e-5)]
        self.assertGreater(max(x)-min(x),1e-6)
        self.assertLess(max(x)-min(x),2e-5)

    def test_16_eta_lambda_order_is_explicit(self):
        for eta in (1e-3,1e-4):
            ref=virtual_probe_limit(5,eta=eta)['Rxy_lambda_to_0_h_over_q2']
            finite=probe_strength_response(5,eta=eta,probe_scale=1e-4)['Rxy_h_over_q2']
            self.assertLess(abs(ref-finite),2e-7)

    def test_17_sign_reversal(self):
        p=virtual_probe_limit(7,mass=+1)
        self.assertAlmostEqual(p['Rxy_lambda_to_0_h_over_q2'],-self.ref['Rxy_lambda_to_0_h_over_q2'],places=10)
        self.assertAlmostEqual(p['two_terminal_T_lr'],self.ref['two_terminal_T_lr'],places=10)

    def test_18_trivial_insulator_fail_closed(self):
        p=virtual_probe_limit(7,mass=3)
        self.assertIsNone(p['Rxy_lambda_to_0_h_over_q2'])
        self.assertFalse(p['limit_defined'])

    def test_19_metalllic_negative_control(self):
        self.assertTrue(self.metal['limit_defined'])
        self.assertLess(abs(self.metal['Rxy_lambda_to_0_h_over_q2']),.3)
        self.assertGreater(abs(self.metal['Rxy_lambda_to_0_h_over_q2']-1),.6)

    def test_20_no_measurement_if_artificially_no_drive_current(self):
        d=virtual_probe_limit(7,min_drive_current=2)
        self.assertIsNone(d['Rxy_lambda_to_0_h_over_q2'])

    def test_21_undetectable_voltage_probe_refused(self):
        d=virtual_probe_limit(7,min_probe_transmission=1000)
        self.assertIsNone(d['Rxy_lambda_to_0_h_over_q2'])
        self.assertIsNone(d['floating_voltages_lambda_to_0'])

    def test_22_width_extreme_is_bounded(self):
        x=virtual_probe_limit(7,drive_width=3,probe_width=2)
        y=virtual_probe_limit(7,drive_width=3,probe_width=7)
        self.assertTrue(x['limit_defined'] and y['limit_defined'])
        self.assertLess(abs(x['Rxy_lambda_to_0_h_over_q2']-y['Rxy_lambda_to_0_h_over_q2']),1e-7)
        self.assertGreater(abs(x['probe_weight_totals'][0]-y['probe_weight_totals'][0]),50)

    def test_23_invalid_geometry(self):
        for kwargs in ({'L':2},{'drive_width':1},{'probe_width':8},{'offsets':(1,)}, {'offsets':(0,10,0,0)}):
            with self.assertRaises(ValueError):virtual_probe_limit(**kwargs)

    def test_24_zero_eta_refused(self):
        with self.assertRaises(ValueError):virtual_probe_limit(eta=0)

    def test_25_nan_refused(self):
        with self.assertRaises(ValueError):virtual_probe_limit(mu=float('nan'))

    def test_26_threshold_positive(self):
        for kwargs in ({'min_probe_transmission':0},{'min_drive_current':-1}):
            with self.assertRaises(ValueError):virtual_probe_limit(**kwargs)

if __name__=='__main__':unittest.main()
