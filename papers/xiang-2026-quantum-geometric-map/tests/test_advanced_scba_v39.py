"""v39 advanced JSON entry must remain scoped to an incomplete Green-sea method."""
import json
import subprocess
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class AdvancedSCBA39(unittest.TestCase):
    def test_unified_json(self):
        run=subprocess.run([sys.executable,str(ROOT/"scripts/reproduce_advanced.py"),
                            "--case","scba-energy-sea-v39","--mesh","17"],
                           cwd=ROOT,text=True,capture_output=True,timeout=80)
        self.assertEqual(run.returncode,0,msg=run.stderr)
        d=json.loads(run.stdout)
        self.assertEqual(d["science_status"],"PARTIAL")
        result=d["result"]
        self.assertLess(result["energy_derivative"]["derivative_vertex_residual"],1e-10)
        self.assertGreater(result["energy_derivative"]["gamma_minus_identity"],.2)
        self.assertIn("full Kubo-Streda Fermi surface",result["open"])
        self.assertGreater(abs(result["metallic"]["derivative_correction"]),.1)
if __name__=="__main__":unittest.main()
