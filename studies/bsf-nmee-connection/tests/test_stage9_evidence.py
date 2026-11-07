"""Protect physical conclusions, including negative findings, from source data."""
import hashlib
import json
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


class Stage9EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = json.loads((ROOT/'results/bilayer_stage9_path.json').read_text())
        cls.coupling = json.loads((ROOT/'results/bilayer_stage9_coupling.json').read_text())

    def test_manifest(self):
        manifest = json.loads((ROOT/'results/STAGE9_MANIFEST.json').read_text())
        for path, digest in manifest['sha256'].items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), digest)

    def test_complete_path_gaps_and_mirror(self):
        for t in (.005,.02):
            rows = [r for r in self.path['path'] if r['t_eV'] == t]
            self.assertEqual(len(rows), 17)
            for a,b in zip(rows, rows[::-1]):
                self.assertGreater(a['sampled_indirect_gap_eV'],.029)
                np.testing.assert_allclose(a['C'], -np.array(b['C']), atol=1e-12)
                self.assertAlmostEqual(a['dipole_over_e_d_per_cell'], -b['dipole_over_e_d_per_cell'], delta=1e-13)
            self.assertLess(rows[0]['C'][0][0]*rows[-1]['C'][0][0],0)
            self.assertLess(rows[0]['dipole_over_e_d_per_cell']*rows[-1]['dipole_over_e_d_per_cell'],0)

    def test_local_ab_expansion_and_refinement(self):
        rows = {r['q']:r for r in self.path['polar_local'] if r['direction'] == [.23,.41]}
        fine = self.path['fine']
        for key in ('dipole_over_e_d_per_cell','trace'):
            coarse = rows[.01][key]+rows[-.01][key]-2*rows[0][key]
            refined = fine['AB_plus'][key]+fine['AB_minus'][key]-2*fine['AB'][key]
            self.assertLess(abs(refined/coarse-1),.001)
            even2 = rows[.02][key]+rows[-.02][key]-2*rows[0][key]
            self.assertAlmostEqual(np.log2(abs(even2/coarse)),2.,delta=.01)

    def test_fixed_stacking_transverse_reversal_and_refinement(self):
        rows = {r['parameters']['tunneling_eV']:r for r in self.coupling['scan']}
        self.assertLess(rows[.005]['C'][1][0],0.)
        self.assertGreater(rows[.02]['C'][1][0],0.)
        for r in self.coupling['scan'][1:]:
            self.assertGreater(r['dipole_over_e_d_per_cell'],0)
            self.assertGreater(r['sampled_indirect_gap_eV'],.031)
        for fine in self.coupling['fine']:
            base = rows[fine['parameters']['tunneling_eV']]
            for key in ('dipole_over_e_d_per_cell','trace','difference'):
                self.assertLess(abs(fine[key]/base[key]-1),.001)

    def test_independent_liouville_agrees_in_numerically_resolved_limit(self):
        data = json.loads((ROOT/'results/bilayer_stage9_limit_check.json').read_text())
        for row in data['models'].values():
            self.assertLess(abs(row['eta_squared_intercept']/row['static']-1),.001)


if __name__ == '__main__':
    unittest.main()
