"""Guard the published scientific conclusions and their numerical uncertainty."""
import hashlib
import json
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((ROOT/'results/bilayer_stage8_audit.json').read_text())
        cls.controls=json.loads((ROOT/'results/bilayer_stage8_controls.json').read_text())

    def test_frozen_evidence_hashes(self):
        manifest=json.loads((ROOT/'results/STAGE8_MANIFEST.json').read_text())
        for path,digest in manifest['sha256'].items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),digest)

    def test_convergence_and_open_gap(self):
        for row in self.data['sliding']:
            self.assertGreater(row['sampled_indirect_gap_eV'],.028)
            self.assertGreater(row['sampled_mu_distance_eV'],0)
        for control in self.data['convergence']:
            base=next(r for r in self.data['sliding'] if r['u']==control['u'])
            for variant in ('fine','half_B'):
                for key in ('dipole_over_e_d_per_cell','trace','difference'):
                    self.assertLess(abs(control[variant][key]/base[key]-1),.001)

    def test_breaking_C3_restores_linear_scalar_terms(self):
        a,b=self.controls['broken_C3_asymptotic'][:2]
        for key in ('dipole_over_e_d_per_cell','trace','difference'):
            self.assertAlmostEqual(np.log2(abs(b[key]/a[key])),1,delta=.02)
        for key in ('dipole_over_e_d_per_cell','trace','difference'):
            self.assertLess(abs(self.controls['broken_C3_fine'][key]/a[key]-1),.001)

    def test_weak_coupling_and_nonadditive_directional_response(self):
        a,b=self.data['channels'][:2]
        self.assertAlmostEqual(b['full']['difference']/a['full']['difference'],4,delta=.05)
        full=a['full']['difference']
        separate=a['same_only']['difference']+a['cross_only']['difference']
        self.assertLess(abs(separate/full),.01)
        self.assertLess(full*a['relative_sign']['difference'],0)
        self.assertLess(np.max(abs(np.array(a['decoupled']['C']))),1e-12)

    def test_flat_spectrum_retains_variation_but_not_Cxx_sign_change(self):
        for t in (.005,.02):
            rows=[r for r in self.data['flat_spectrum'] if r['t_eV']==t]
            self.assertLess(rows[0]['physical']['C'][0][0],0)
            self.assertGreater(rows[-1]['physical']['C'][0][0],0)
            self.assertTrue(all(r['flat']['C'][0][0]<0 for r in rows))
            self.assertGreater(abs(rows[-1]['flat']['C'][0][0]/rows[0]['flat']['C'][0][0]),10)
        for control in self.controls['flat_convergence']:
            base=next(r['flat'] for r in self.data['flat_spectrum']
                      if r['t_eV']==control['t_eV'] and r['u']==control['u'])
            for variant in ('fine','half_B'):
                for key in ('trace','difference'):
                    self.assertLess(abs(control[variant][key]/base[key]-1),.0001)


if __name__=='__main__':
    unittest.main()
