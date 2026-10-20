"""Fail release if archived stage-6 scientific inputs/figures differ."""
import unittest,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class Stage6ReleaseGate(unittest.TestCase):
    def test_frozen_sources_and_numbers(self):
        m=json.loads((ROOT/'results/STAGE6_MANIFEST.json').read_text())
        self.assertEqual(m['scientific_status'],'synthetic_8band_static_only')
        for path,digest in m['sha256'].items():
            actual=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            self.assertEqual(actual,digest,path)
        r=json.loads((ROOT/'results/sliding_stage6_symmetry.json').read_text())
        self.assertLess(abs(r['rotated_shifts']['60']['trace_ratio']+1),.002)
        self.assertLess(r['rotated_shifts']['120']['matrix_residual'],2e-10)
        self.assertTrue(7.3<r['empirical_scale_ratios_u01_over_u005']['trace']<8.5)
        self.assertTrue(1.8<r['empirical_scale_ratios_u01_over_u005']['mixed_nematic']<2.2)
if __name__=='__main__':unittest.main()
