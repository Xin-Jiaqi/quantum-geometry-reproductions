"""The documented disorder/Kubo/Ward CLI cases emit valid, scoped JSON."""
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

class DisorderCliContract(unittest.TestCase):
    def test_new_cases_json(self):
        for case in ("scba-ra","kubo-spectral","finite-qomega-ward"):
            with self.subTest(case=case):
                run=subprocess.run(
                    [sys.executable,str(ROOT/"scripts"/"reproduce_advanced.py"),
                     "--case",case,"--mesh","12"],
                    cwd=ROOT,text=True,capture_output=True,check=False)
                self.assertEqual(run.returncode,0,msg=run.stderr)
                payload=json.loads(run.stdout)
                self.assertEqual(payload["case"],case)
                self.assertEqual(payload["science_status"],"PARTIAL")
                self.assertIn("result",payload)
if __name__=="__main__":unittest.main()
