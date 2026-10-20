"""Unified advanced science-status contract for v38 Hall-sea benchmark."""
import json
import subprocess
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
class StredaAdvancedCase(unittest.TestCase):
    def test_case_stays_clean_and_partial(self):
        p=subprocess.run([sys.executable,str(ROOT/"scripts/reproduce_advanced.py"),
            "--case","clean-streda-sea-v38"],cwd=ROOT,text=True,capture_output=True,timeout=45)
        self.assertEqual(p.returncode,0,msg=p.stderr)
        payload=json.loads(p.stdout)
        self.assertEqual(payload['science_status'],'PARTIAL')
        d=payload['result']
        self.assertEqual(d['scope'],'clean_antisymmetric_sea_only')
        self.assertAlmostEqual(d['fhs_chern'],-1.,delta=1e-12)
        self.assertAlmostEqual(d['sea_sigma_over_e2h'],-.9750936556674081,delta=1e-10)
        self.assertIn('complete_metallic_disorder',d['not_implemented'])
if __name__=="__main__":unittest.main()
