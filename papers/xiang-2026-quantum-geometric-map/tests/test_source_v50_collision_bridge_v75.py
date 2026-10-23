"""v75: optional full-repository original-source v49/v50 bridge.

Standalone handoff intentionally omits old copyrighted source and its
pre-v60 code; these tests are skipped there, not reported as passed.
"""
import unittest
from dataclasses import replace
try:
    from qgmap.sm_literal_spin_phe_v49 import SourceConfig, model
    from qgmap.sm_two_color_liouville_v50 import Tone, two_color
    HAS_FULL_SOURCE=True
except ImportError:
    HAS_FULL_SOURCE=False
from qgmap.selected_electric_pole_v74 import Spec, literal_eq20, printed_eq31
from qgmap.elastic_spin_collision_v75 import Spec as CollisionSpec, steady

@unittest.skipUnless(HAS_FULL_SOURCE,'standalone v60+ snapshot lacks original v49/v50 source modules')
class SourceV50VsElasticV75(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=SourceConfig(nr=44,nphi=64,temperature=.03)
        cls.s=Spec(nr=44,nphi=64,temperature=.03)
        cls.lit=model(cls.c,'main')['xy']

    def test_original_sm_literal_eq20_independent_v74(self):
        self.assertLess(abs(self.lit['SM_Eq20_current_physical_B']-literal_eq20(self.s)),3e-11)

    def test_original_sm_printed_eq31_relative_sign_not_fixed(self):
        self.assertLess(abs(self.lit['SM_printed_Eq31_current_physical_B']-printed_eq31(self.s)),3e-11)
        self.assertGreater(abs(self.lit['SM_Eq20_current_physical_B']-self.lit['SM_printed_Eq31_current_physical_B']),.004)

    def test_original_v50_selected_pole_matches_eq20_not_eq31(self):
        a=two_color(self.c,Tone(eta_a=.001,eta_b=.001,dk=2e-5),'main','xy')
        self.assertLess(abs(a['eta_a_times_total'].real-self.lit['SM_Eq20_current_physical_B']),2e-6)
        self.assertGreater(abs(a['eta_a_times_total'].real-self.lit['SM_printed_Eq31_current_physical_B']),.004)

    def test_dirac_scattering_transport_tau_not_quantum_tau(self):
        tau=steady(CollisionSpec(beta=.125))['tau_current'].real
        self.assertGreater(abs(tau-1/.07),10.)
        self.assertLess(abs(tau-2/.07),1.)

if __name__=='__main__':unittest.main()
