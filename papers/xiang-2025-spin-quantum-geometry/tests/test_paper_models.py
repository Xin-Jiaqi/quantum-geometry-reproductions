"""Offline paper-model acceptance without publisher source files."""
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import paper_models as m

class PaperModelTests(unittest.TestCase):
    def test_dirac_spectrum(self):
        self.assertAlmostEqual(np.diff(m.tilted_bands(0.))[0],.08,places=12)
    def test_tilted_analytic_spin_qgt(self):
        for k in ((.01,.02),(-.07,.04),(.0,.0)):
            g,_=m.spin_qgt(*k)
            self.assertAlmostEqual(g,m.sm_metric(*k),places=10)
    def test_sm_metric_exact_angular(self):
        for w in (20.,40.,60.):
            self.assertAlmostEqual(m.sm_metric_independent_2d(w),m.sm_metric_from_eq55(w),delta=3e-11)
    def test_sm_source_discrepancy_is_retained(self):
        self.assertGreater(m.sm_metric_from_eq55(40.)-m.sm_metric_printed(40.),.02)
    def test_curvature_origin(self):
        self.assertAlmostEqual(m.sm_curvature(0,0),.2/(4*.04**2),places=9)
    def test_fig1a_two_valleys_and_no_drude_in_gap(self):
        p=.6582119569*100.;alpha=.1*p
        self.assertEqual(m.fig1a_ndsc(0,p,p,alpha,10.),0.)
        self.assertNotEqual(m.fig1a_ndsc(20.,p,p,alpha,10.),0.)
    def test_fig1a_spin_metric_matrix(self):
        for valley in (-1,1):
            for kx,ky in ((.12,-.07),(-.31,.26)):
                a=m.fig1a_spin_metric(kx,ky,valley,65.8,65.8,6.58,10.)
                b=m.fig1a_matrix_metric(kx,ky,valley,65.8,65.8,6.58,10.)
                self.assertAlmostEqual(a,b,delta=1e-11)
    def test_fig1b_branches_not_merged(self):
        r=m.main_fig1b(0.)
        self.assertLess(r['intrinsic_maintext_natural'],0.)
        self.assertGreater(r['intrinsic_literal_SM_integrand_natural'],0.)
        self.assertNotAlmostEqual(r['intrinsic_maintext_natural'],r['intrinsic_SM70_natural'],delta=1e-10)
    def test_fig1b_drude_in_gap(self):
        self.assertEqual(m.main_fig1b(0.)['drude_SM67_natural'],0.)
    def test_frequency_scope_guard(self):
        for w in (0.,80.,100.):
            with self.assertRaises(ValueError):m.sm_metric_printed(w)
    def test_report_status(self):
        d=m.reproduce()
        self.assertEqual(d['status'],'PARTIAL')
        self.assertEqual(d['scientific_gates']['metal_DC_protocol'],'NOT_VALIDATED')
        self.assertNotAlmostEqual(d['sm_fig1b'][0]['printed_eq56'],d['sm_fig1b'][0]['integrated_eq55'],delta=.001)

if __name__=='__main__':unittest.main()
