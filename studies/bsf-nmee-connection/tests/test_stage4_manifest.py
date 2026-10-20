"""Frozen stage-4 outputs, source integrity and scientific limits."""
import json,hashlib,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class TestStage4Manifest(unittest.TestCase):
    def test_manifest(self):
        manifest=json.loads((ROOT/'results/STAGE4_MANIFEST.json').read_text())
        for rel,digest in manifest['sha256'].items():
            p=ROOT/rel
            self.assertTrue(p.is_file(),rel)
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),digest,rel)
    def test_frozen_numerical_physics(self):
        a=json.loads((ROOT/'results/sliding_stage4_mechanism.json').read_text())
        p=a['hopping_channel_competition']['0.1']
        q=a['hopping_channel_competition']['0.4']
        self.assertLess(p['full']['gamma'],0)
        self.assertGreater(p['cross_only']['gamma'],0)
        self.assertLess(p['inclusion_exclusion_nonadditivity'],0)
        self.assertGreater(q['inclusion_exclusion_nonadditivity'],0)
        self.assertGreater(a['relative_channel_sign_controls']['0.1']['cross_minus_gamma'],0)
        self.assertGreater(a['parameter_robustness']['baseline']['min_sampled_mu_distance_eV'],.001)
        self.assertGreater(a['parameter_robustness']['weak_hopping_10meV']['min_sampled_mu_distance_eV'],.005)
        self.assertGreater(a['parameter_robustness']['weak_hopping_10meV']['zero_midpoint'],a['parameter_robustness']['baseline']['zero_midpoint'])
        for x in a['independent_static_kernel_and_quadrature_checks'].values():
            self.assertLess(x['mesh_change'],5e-9)
if __name__=='__main__':unittest.main()
