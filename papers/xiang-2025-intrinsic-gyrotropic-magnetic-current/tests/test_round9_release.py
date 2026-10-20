"""Metadata and packaging contracts; not a claim that Fig.2(d) matches."""
import json
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from verify_round9 import preflight

class TestReleaseCandidate(unittest.TestCase):
    def test_documented_partial_scope(self):
        contract=preflight()
        self.assertEqual(contract['scientific_status'],'PARTIAL')
        self.assertFalse(contract['allow_empirical_multiplier'])

    def test_expected_two_response_protocols(self):
        c=preflight()
        self.assertEqual(c['response_branches'],[
            'paper-defined-low-frequency', 'independent-full-frequency'])

    def test_author_data_blocked(self):
        c=preflight()
        self.assertEqual(c['scope']['original_fig2d_absolute'],'BLOCKED_AUTHOR_DATA')
        self.assertEqual(c['author_contact'],'NOT_SENT')

    def test_magnetic_toolkit_license_retained(self):
        text=(ROOT/'data/TOOLKIT_BSD3_LICENSE.txt').read_text()
        self.assertIn('BSD',text)
        self.assertTrue((ROOT/'data/toolkit_mpg_generators_v018.txt').is_file())

    def test_prior_frozen_has_required_checkpoints(self):
        c=preflight()
        f=json.loads((ROOT/c['frozen_baseline_round8']).read_text())
        self.assertIn('prior_frozen', f)
        self.assertIn('expected',f)

if __name__=='__main__':unittest.main()
