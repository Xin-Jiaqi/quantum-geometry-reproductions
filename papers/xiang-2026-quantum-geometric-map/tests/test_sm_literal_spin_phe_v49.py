"""Source-certified v49 literal SM Eq20/31/32 vs main/SM Dirac helicity tests."""
import unittest
import numpy as np
from qgmap.sm_literal_spin_phe_v49 import SourceConfig,model,compare_source,geometry,SX,SY,ID

class SourceCertifiedPHEV49(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=SourceConfig()
        cls.a=compare_source(cls.s)

    def test_literal_SM_20_vs_SM_31_sign_both_helicities(self):
        for t in (self.a['main'],self.a['supplementary']):
            for ab in ('xy','yx'):
                self.assertLess(abs(t[ab]['Eq20_plus_Eq31_residual']),1e-12)
                self.assertGreater(abs(t[ab]['SM_printed_Eq31_current_physical_B']),1e-4)

    def test_explicit_symmetry_of_printed_31_and_32(self):
        for t in (self.a['main'],self.a['supplementary']):
            q=.5*(t['xy']['SM_printed_Eq31_current_physical_B']+t['yx']['SM_printed_Eq31_current_physical_B'])
            self.assertAlmostEqual(q,t['SM_printed_Eq32_symmetric'],delta=1e-14)
            self.assertLess(abs(t['xy']['SM_printed_Eq31_current_physical_B']-t['yx']['SM_printed_Eq31_current_physical_B']),1e-8)

    def test_main_and_SM_Hamiltonians_are_opposite_fixed_spin_helicities(self):
        x,y,b=.27,-.43,.125
        H_main=b*x*ID+x*SY-y*SX
        H_sm=b*x*ID+y*SX-x*SY
        sz=np.diag([1.,-1.])
        np.testing.assert_allclose(sz@H_main@sz,H_sm,atol=1e-14)
        np.testing.assert_allclose(sz@SX@sz,-SX,atol=1e-14)
        self.assertGreater(np.max(abs(H_main-H_sm)),.5)

    def test_main_SM_sign_flip_in_direct_Z_geometry(self):
        a=self.a['main']['SM_printed_Eq32_symmetric']
        b=self.a['supplementary']['SM_printed_Eq32_symmetric']
        self.assertAlmostEqual(a,-b,delta=1e-13)
        self.assertGreater(a,0)
        self.assertLess(b,0)

    def test_literal_SM_62_and_63_match_with_SM_helicity(self):
        t=self.a['supplementary']
        self.assertLess(abs(t['Eq32_to_Eq11_ratio']-1),2e-6)
        self.assertGreater(t['xy']['SM_Eq20_current_physical_B'],0)  # Eq20 NOT Eq31
        self.assertLess(t['SM_printed_Eq32_symmetric'],0)

    def test_published_main_Eq10_vs_Eq11_have_opposite_formula_sign(self):
        t=self.a['main']
        self.assertLess(abs(t['Eq32_to_Eq11_ratio']+1),2e-6)
        self.assertGreater(t['SM_printed_Eq32_symmetric'],0)
        self.assertLess(t['main_Eq11_geometry_tanh_extension'],0)

    def test_independent_analytic_main_Eq11_identity(self):
        for beta in (.025,.125,.3,.6):
            exact=(beta*beta+2*np.sqrt(1-beta*beta)-2)/(4*np.pi*beta**3)
            self.assertAlmostEqual(geometry(beta),exact,delta=2e-9)

    def test_beta_scans_are_neither_special_beta_nor_sign_fit(self):
        for beta in (.025,.08,.3,.6):
            s=SourceConfig(beta=beta,nr=96,nphi=128)
            main=model(s,'main')
            sm=model(s,'sm')
            self.assertLess(abs(main['Eq32_to_Eq11_ratio']+1),2e-5)
            self.assertLess(abs(sm['Eq32_to_Eq11_ratio']-1),2e-5)

    def test_gauge_independence_even_under_different_band_phases(self):
        ref=self.a['main']
        twisted=model(self.s,'main',(0.23,-.81))
        for c in ('xy','yx'):
            for key in ('SM_Eq20_current_physical_B','SM_printed_Eq31_current_physical_B'):
                self.assertAlmostEqual(twisted[c][key],ref[c][key],delta=1e-12)

    def test_electron_hole_parity_low_temperature_model(self):
        pos=self.a['supplementary']
        neg=model(SourceConfig(mu=-.3),'sm')
        self.assertAlmostEqual(pos['SM_printed_Eq32_symmetric'],-neg['SM_printed_Eq32_symmetric'],delta=1e-9)

    def test_adiabatic_pole_partial_fraction_identity(self):
        # Literal SM Eq15 rational identity BEFORE one substitutes 1/eta_a by tau.
        for delta in (.31,1.27):
            for a,b in ((.0002,.0003),(.02,.11)):
                A=-delta+1j*(a+b)
                B=-delta+1j*b
                lhs=1/(A*1j*a)
                rhs=(1/(1j*a)-1/A)/B
                self.assertAlmostEqual(lhs.real,rhs.real,delta=1e-8)
                self.assertAlmostEqual(lhs.imag,rhs.imag,delta=1e-8)

    def test_fail_closed(self):
        with self.assertRaises(ValueError):SourceConfig(beta=1.)
        with self.assertRaises(ValueError):SourceConfig(temperature=-1.)
        with self.assertRaises(ValueError):SourceConfig(nr=5)
        with self.assertRaises(ValueError):model(helicity='unknown')
        with self.assertRaises(ValueError):geometry(-.125)

if __name__=='__main__':unittest.main()
