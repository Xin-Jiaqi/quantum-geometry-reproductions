"""v37 API, exact provenance and published-figure model-boundary regression."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/reproduce_advanced.py"

def call(case,mesh=None):
    cmd=[sys.executable,str(SCRIPT),"--case",case]
    if mesh is not None:cmd.extend(["--mesh",str(mesh)])
    completed=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=75)
    if completed.returncode:raise AssertionError(completed.stderr)
    payload=json.loads(completed.stdout)
    assert payload["science_status"]=="PARTIAL"
    return payload["result"]

class V37PublishedModel(unittest.TestCase):
    def test_continuum_kernel(self):
        r=call("continuum-dirac-v24")
        self.assertLess(r["projector_eigensolver_abs_diff"],1e-10)
        self.assertLess(r["sm20_angular_abs_diff"],1e-8)
        self.assertAlmostEqual(r["sm20_direct_per_bbar_tau"]["real"],0.0024899354131046854,delta=1e-8)

    def test_si_two_distinct_protocol_branches(self):
        r=call("fig2-si-v25")
        self.assertAlmostEqual(r["literal_hbar"]["minus_1e4_sigma_over_e2_hbar"],11.02083611497046,delta=1e-8)
        self.assertAlmostEqual(r["single_2pi_diagnostic"]["minus_1e4_sigma_over_e2_hbar"],1.754020544703229,delta=1e-8)
        self.assertAlmostEqual(r["ratio"],2*math.pi,delta=1e-10)

    def test_independent_joule_liouville(self):
        r=call("si-liouville-v26")
        self.assertLess(r["relative_magnitude_error"],5e-5)
        self.assertAlmostEqual(r["literal_eq11_magnitude_S"],2.6823765071653156e-7,delta=3e-11)

    def test_digitized_curve_is_approximate(self):
        r=call("fig2-curve-v27")
        for label,T in (("T20",20),("T40",40)):
            self.assertAlmostEqual(r["fits"][label]["effective_T_K"],T,delta=.7)
            self.assertAlmostEqual(r["fits"][label]["amplitude"],1.75,delta=.012)
        self.assertAlmostEqual(r["platform_mean"],1.7495,delta=.012)
        self.assertLess(abs(r["ratio_to_literal"]-1/(2*math.pi)),.003)
        self.assertIn("axis_calibration_sensitivity",r)

    def test_snapshot_provenance_exact_sha256_and_git_blobs(self):
        reference=json.loads((ROOT/"results/frozen/v37_manifest.json").read_text())
        self.assertEqual(len(reference),13)
        for name,expected in reference.items():
            relative=(
                "tests/" if name.startswith("test_") else
                "src/qgmap/" if name.endswith(".py") else
                "results/raw/"
            )
            data=(ROOT/relative/name).read_bytes()
            with self.subTest(name=name):
                self.assertEqual(hashlib.sha256(data).hexdigest(),expected["sha256"])
                hdr=f"blob {len(data)}\0".encode()
                self.assertEqual(hashlib.sha1(hdr+data).hexdigest(),expected["git_blob_sha1"])

    def test_invalid_mesh_fails_closed(self):
        for case,n in (("continuum-dirac-v24",10),("si-liouville-v26",10),("fig2-si-v25",60),("fig2-curve-v27",60)):
            c=subprocess.run([sys.executable,str(SCRIPT),"--case",case,"--mesh",str(n)],cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(c.returncode,0,msg=case)

if __name__=="__main__":unittest.main()
