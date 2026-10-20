"""Frozen model numbers from independently verified V27 baseline.

No normalization is tuned to printed curves. Updating references requires an
explicit scientific review rather than running a snapshot overwrite command.
"""
import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import paper_models as m

class FrozenRelease(unittest.TestCase):
    def test_numerical_and_scientific_status(self):
        ref=json.loads((ROOT/'results/REFERENCE_V28.json').read_text())
        d=m.reproduce()
        mu0=m.main_fig1b(0.)
        got={
            'main_fig1a_mu20_ndsc_natural':d['fig1a']['ndsc_mu20'],
            'main_fig1b_mu0_maintext_natural':mu0['intrinsic_maintext_natural'],
            'main_fig1b_mu0_sm70_natural':mu0['intrinsic_SM70_natural'],
            'main_fig1b_mu0_static_source_natural':mu0['independent_density_matrix_static_natural'],
            'sm_fig1b_omega40_printed_eq56':m.sm_metric_printed(40.),
            'sm_fig1b_omega40_integrated_eq55':m.sm_metric_independent_2d(40.),
            'sm_fig1b_omega40_berry_eq58':m.sm_berry_printed(40.),
        }
        self.assertEqual(ref['scientific_status'],'PARTIAL')
        self.assertEqual(d['status'],ref['scientific_status'])
        for key,actual in got.items():
            with self.subTest(key=key):
                target=ref[key]
                self.assertAlmostEqual(actual,target,delta=1e-10*max(abs(target),1e-6))

if __name__=='__main__':unittest.main()
