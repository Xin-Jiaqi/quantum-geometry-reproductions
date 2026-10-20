"""Verify submitted stage-3 source/data integrity, not physical validation."""
import unittest, json, hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
class TestStage3Manifest(unittest.TestCase):
    def test_manifest_matches_all_released_sources(self):
        manifest=json.loads((ROOT/'results/STAGE3_MANIFEST.json').read_text())
        self.assertEqual(manifest['files'],len(manifest['sha256']))
        for relative,digest in manifest['sha256'].items():
            p=ROOT/relative
            self.assertTrue(p.is_file(),relative)
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),digest,relative)
    def test_frozen_sample_and_root_are_consistent(self):
        d=json.loads((ROOT/'results/sliding_stage3_sign_reversal.json').read_text())
        scan={x['path_parameter']:x for x in d['scan']}
        self.assertLess(scan[0.1]['gamma'],0)
        self.assertGreater(scan[0.2]['gamma'],0)
        self.assertLess(scan[0.1]['sampled_direct_gap_eV'],scan[0.4]['sampled_direct_gap_eV'])
        self.assertGreater(min(x['sampled_min_mu_to_band_eV'] for x in d['scan']),.001)
        lo,hi=d['zero_crossing']['root_bracket'];self.assertLess(lo,.198)
        self.assertGreater(hi,.197)
        for val in d['mesh_convergence'].values():self.assertLess(val['absolute_change'],3e-9)
        for val in d['mirror_oddness'].values():self.assertLess(abs(val['odd_residual']),1e-11)
if __name__=='__main__':unittest.main()
