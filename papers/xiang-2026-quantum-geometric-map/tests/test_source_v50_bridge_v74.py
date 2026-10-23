"""v74 independent authentic-SM v50 finite-difference cross-oracle.

Requires the repository's original v49/v50/Dirac integration modules;
portable v60+ handoff correctly marks these source-archive tests SKIPPED.
Full GitHub Actions runs them (the repository contains v49/v50).
"""
import unittest
from qgmap.selected_electric_pole_v74 import Spec,channels,literal_eq20
try:
    from qgmap.sm_two_color_liouville_v50 import Tone,two_color
    from qgmap.sm_literal_spin_phe_v49 import SourceConfig,model
    HAS_ORIGINAL_SOURCE=True
except ImportError:
    HAS_ORIGINAL_SOURCE=False

@unittest.skipUnless(HAS_ORIGINAL_SOURCE,'v49/v50 original source modules live in GitHub complete repo, not standalone v60+ handoff')
class V74OriginalV50IndependentOracle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(nr=64,nphi=88,eta_e=.003,eta_b=.001)
        cls.c=SourceConfig(nr=64,nphi=88,beta=.125,mu=.3,temperature=.03,cutoff=2.5)
        cls.t=Tone(eta_a=.003,eta_b=.001,dk=1e-5)
        cls.old=two_color(cls.c,cls.t,helicity='main',direction='xy')
        cls.now=channels(cls.s)

    def test_v50_full_finite_difference_matches_new_frechet(self):
        self.assertLess(abs(self.old['eta_a_times_total']-self.now['eta_e_scaled_full']),4e-7)

    def test_v50_B_then_E_independent_k_fd(self):
        self.assertLess(abs(self.old['eta_a_times_b_then_e']-self.now['eta_e_scaled_b_then_e']),4e-7)

    def test_v50_E_then_B_independent_k_fd(self):
        self.assertLess(abs(self.old['eta_a_times_e_then_b']-self.now['eta_e_scaled_e_then_b']),4e-7)

    def test_v50_selected_explicit_diagonal_pole(self):
        old=two_color(self.c,self.t,pole_only=True)['eta_a_times_physical_B_current_residue']
        self.assertLess(abs(old-self.now['eta_e_scaled_intraband_e_pole']),4e-7)

    def test_v49_literal_published_SM_Eq20_branch(self):
        old=model(self.c,'main')['xy']['SM_Eq20_current_physical_B']
        self.assertLess(abs(old-literal_eq20(self.s)),5e-7)
        self.assertLess(abs(self.now['eta_e_scaled_full']-old),1e-6)

    def test_v49_printed_Eq31_relative_sign_remains_unfixed(self):
        old=model(self.c,'main')['xy']['SM_printed_Eq31_current_physical_B']
        self.assertLess(abs(old+literal_eq20(self.s)),5e-7)
        self.assertGreater(abs(self.now['eta_e_scaled_full']-old),.004)

if __name__=='__main__':unittest.main()
