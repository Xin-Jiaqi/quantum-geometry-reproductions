"""v35 published method interfaces against frozen v20 independent model numbers."""
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/"results"/"frozen"/"v35_real_frequency_reference.json").read_text())

def run_case(case,mesh=None):
    args=[sys.executable,str(ROOT/"scripts"/"reproduce_advanced.py"),"--case",case]
    if mesh is not None:args.extend(["--mesh",str(mesh)])
    completed=subprocess.run(args,cwd=ROOT,text=True,capture_output=True,timeout=90)
    if completed.returncode: raise AssertionError(completed.stderr)
    return json.loads(completed.stdout)

class V35Contracts(unittest.TestCase):
    def test_exact_json_and_reference(self):
        d=run_case("exact-kubo-v20")
        self.assertEqual(d["science_status"],"PARTIAL")
        r=d["result"]
        self.assertEqual(r["L"],5)
        self.assertEqual(r["hilbert_space_dim"],50)
        self.assertLess(abs(r["thermodynamic"]["grand_canonical_compressibility"]-DATA["exact"]["compressibility"]),1e-10)
        self.assertLess(abs(r["directions"]["y"]["missing_contact_abs"]-DATA["exact"]["missing_contact_y_abs"]),1e-9)
        for x in ("x","y"):
            self.assertLess(r["directions"][x]["ward_density_abs"],1e-10)
            self.assertLess(r["directions"][x]["ward_current_abs"],1e-10)
    def test_scba_limit_json_and_frozen(self):
        d=run_case("limit-order-v20")
        self.assertEqual(d["science_status"],"PARTIAL")
        seq=d["result"]["values"]
        self.assertEqual(len(seq),3)
        for idx,x in enumerate(seq):
            self.assertAlmostEqual(x["eta"],DATA["scba"]["eta"][idx],places=12)
            self.assertLess(abs(x["fixed_probe_regulator_leak_abs"]-DATA["scba"]["fixed_probe_leak"][idx]),3e-8)
            self.assertLess(abs(x["uniform_static_regulator_leak_abs"]-DATA["scba"]["uniform_static_leak"][idx]),3e-8)
            self.assertLess(x["fixed_probe_ward_abs"],2e-9)
        self.assertLess(seq[-1]["fixed_probe_regulator_leak_abs"],seq[0]["fixed_probe_regulator_leak_abs"])
    def test_different_notions_of_mesh_rejected(self):
        d=run_case("exact-kubo-v20",4)
        self.assertEqual(d["result"]["L"],4)
        reject=subprocess.run([sys.executable,str(ROOT/"scripts"/"reproduce_advanced.py"),
           "--case","exact-kubo-v20","--mesh","14"],cwd=ROOT,capture_output=True,text=True)
        self.assertNotEqual(reject.returncode,0)

if __name__=="__main__":unittest.main()
