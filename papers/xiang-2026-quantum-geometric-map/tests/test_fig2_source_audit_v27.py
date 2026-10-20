"""Standalone tests for derived published-curve shape, not author raw datasets."""
import unittest,sys,math,json
from pathlib import Path
import numpy as np
from scipy import constants as c
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from qgmap.fig2_source_audit_v27 import load_points,fit_branch,thermal_shape,px_to_values,ratio_to_literal,calibration_sensitivity
from qgmap.physical_si_liouville_v26 import SIConfig,eq20_physical,eq11_SI_magnitude
DATA=ROOT/'results/raw/v27_fig2b_digitized.csv'
RAW=ROOT/'results/raw/v27_curve_shape_and_amplitude.json'
class TestPublishedCurveAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pts=load_points(DATA)
        cls.out=json.loads(RAW.read_text())
    def test_curve_samples_present(self):
        for name in ('T20','T40'):self.assertGreaterEqual(len(self.pts[name]),25)
    def test_sorted_positive_branch(self):
        for arr in self.pts.values():
            a=np.array(arr);self.assertTrue(np.all(np.diff(a[:,0])>0));self.assertTrue(np.all(a[:,0]>0))
    def test_20K_temperature(self):self.assertLess(abs(fit_branch(self.pts['T20'])['effective_T_K']-20),1.2)
    def test_40K_temperature(self):self.assertLess(abs(fit_branch(self.pts['T40'])['effective_T_K']-40),1.5)
    def test_20K_fit_residual(self):self.assertLess(fit_branch(self.pts['T20'])['rms_ordinate'],.01)
    def test_40K_fit_residual(self):self.assertLess(fit_branch(self.pts['T40'])['rms_ordinate'],.01)
    def test_plateau_height_20(self):self.assertLess(abs(fit_branch(self.pts['T20'])['amplitude']-1.75),.05)
    def test_plateau_height_40(self):self.assertLess(abs(fit_branch(self.pts['T40'])['amplitude']-1.75),.05)
    def test_ratio_compares_to_two_pi_not_exact(self):
        fitA=np.mean([fit_branch(x)['amplitude'] for x in self.pts.values()]);ratio,_=ratio_to_literal(fitA)
        self.assertLess(abs(ratio*2*math.pi-1),.012)
        self.assertGreater(abs(ratio-1),.5)
    def test_calibration_envelope_contains_printed_temperature(self):
        r=calibration_sensitivity(DATA)
        for t,name in [(20,'T20'),(40,'T40')]:
            self.assertLessEqual(r[name]['T_K_range'][0],t)
            self.assertGreaterEqual(r[name]['T_K_range'][1],t)
    def test_flat_mu_zero_and_symmetry(self):
        for T in [20,40]:
            self.assertEqual(thermal_shape([0],T)[0],0)
            self.assertLess(np.max(abs(thermal_shape(np.array([-3,-1,1,3]),T)+thermal_shape(np.array([3,1,-1,-3]),T))),1e-14)
    def test_source_does_not_claim_raw_publisher_data(self):
        self.assertIn('not authors raw data',self.out['curve_extraction'])
        self.assertTrue('PARTIAL' in self.out['status'])
    def test_frozen_values(self):
        for name in ('T20','T40'):
            expected=self.out['fit'][name]; actual=fit_branch(self.pts[name]);
            self.assertAlmostEqual(expected['amplitude'],actual['amplitude'],places=9)
            self.assertAlmostEqual(expected['effective_T_K'],actual['effective_T_K'],places=7)
    def test_independent_SI_bands_finite_T(self):
        for T,mu in [(20,2),(40,5)]:
            p=SIConfig(mu_eV=mu/1000,thermal_eV=c.k*T/c.e,nr=128,nphi=128)
            d=abs(eq20_physical(p).real/p.E_Vm)
            self.assertLess(abs(d/eq11_SI_magnitude(p)-1),.0005)
    def test_analytic_frozen_SI_thermal_batch(self):
        self.assertLess(max(t['relative_deviation'] for t in self.out['independent_finiteT_eq20_check']),.001)
if __name__=='__main__':unittest.main()
