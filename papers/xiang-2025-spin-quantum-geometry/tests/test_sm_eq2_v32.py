"""Nondegenerate Eq.2 four source regression, migrated from offline V14."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from methods.sm_eq2_v32 import (Model,probe,eq2_four_sources,eta0_eq2_and_printed_eq6,
    rho1_band,printed_eq6_triple_independent)
class SMFourSourceTests(unittest.TestCase):
    def test_01_two_orbital(self):self.assertLess(probe(2)['eq2_vs_orbital_max'],2e-8)
    def test_02_three_orbital(self):self.assertLess(probe(3)['eq2_vs_orbital_max'],2e-8)
    def test_03_four_orbital(self):self.assertLess(probe(4)['eq2_vs_orbital_max'],2e-8)
    def test_04_five_orbital(self):self.assertLess(probe(5)['eq2_vs_orbital_max'],2e-8)
    def test_05_three_band_triple(self):self.assertGreater(probe(3)['triple_max'],5e-4)
    def test_06_four_band_triple(self):self.assertGreater(probe(4)['triple_max'],5e-3)
    def test_07_three_band_missing_term(self):self.assertLess(probe(3)['eq2_minus_eq6_max'],1e-13)
    def test_08_five_band_missing_term(self):self.assertLess(probe(5)['eq2_minus_eq6_max'],1e-13)
    def _gauge(self,n):
        a=probe(n);b=probe(n,phi_slope=tuple(.23*(i+1) for i in range(n)),
                             phi_offset=tuple(.37*i for i in range(n)))
        for key in ('delta_spin_trace','eq2_offdiag_spin_trace',
                    'eq6_offdiag_spin_trace','full_orbital_spin_trace'):
            self.assertAlmostEqual(a[key],b[key],delta=2e-9)
        self.assertLess(b['eq2_vs_orbital_max'],2e-8)
    def test_09_gauge3(self):self._gauge(3)
    def test_10_gauge4(self):self._gauge(4)
    def test_11_gauge5(self):self._gauge(5)
    def test_12_eta_asymptotic(self):
        m=Model(4);c=eta0_eq2_and_printed_eq6(m,.23)
        mask=~np.eye(4,dtype=bool)
        errors=[]
        for eta in (.08,.04,.02,.01):
            rho=eq2_four_sources(m,.23,eta).sum(0)
            errors.append(float(np.max(abs((rho-c['eq4_leading']/eta-c['eq2_C'])[mask]))))
        self.assertTrue(all(a>b for a,b in zip(errors,errors[1:])))
        self.assertLess(errors[-1],1.5e-4)
    def test_13_flat_occ_no_extra(self):
        m=Model(4);m.slopes[:]=0.;m.curvature[:]=0.
        x=eta0_eq2_and_printed_eq6(m,.23)
        self.assertLess(np.max(abs(x['extra'])),1e-14)
    def test_14_second_seed(self):
        x=probe(N=4,seed=80888,k=-.31,eta=.015,mu=-.15,T=.31)
        self.assertLess(x['eq2_vs_orbital_max'],2e-8)
        self.assertGreater(x['triple_max'],1e-5)
    def test_15_rho1_hermitian(self):
        d=Model(5).data(.43,mu=.22,T=.3);x=rho1_band(d,.017)
        self.assertLess(np.max(abs(x-x.conj().T)),1e-13)
    def test_16_source_sum_hermitian(self):
        x=eq2_four_sources(Model(4),.23,.01).sum(0)
        self.assertLess(np.max(abs(x-x.conj().T)),1e-11)
    def test_17_nondegenerate_energy_gaps(self):
        for n in (3,4,5):self.assertGreater(np.min(np.diff(Model(n).data(.23)['en'])),.5)
    def test_18_step_stability(self):
        m=Model(4);a=eta0_eq2_and_printed_eq6(m,.23,h=2e-5)
        b=eta0_eq2_and_printed_eq6(m,.23,h=1e-5)
        self.assertLess(np.max(abs(a['eq2_C']-b['eq2_C'])),3e-9)
    def test_19_independent_triple_formula(self):
        for n in (2,3,4,5):
            m=Model(n)
            a=eta0_eq2_and_printed_eq6(m,.23)['triple']
            b=printed_eq6_triple_independent(m,.23)
            self.assertLess(np.max(abs(a-b)),1e-14)
