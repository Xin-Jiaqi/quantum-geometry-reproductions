"""Stage8p: require a NONZERO sign-correct 3-route static-response benchmark."""
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/"hpc-stage8i"
sys.path.insert(0,str(ROOT))
spec=importlib.util.spec_from_file_location("stage8p_nonzero_benchmark",ROOT/"stage8p_nonzero_benchmark.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class NonzeroSourceResponseTests(unittest.TestCase):
    def test_nonzero_analytic_resolvent_projector_wilson_agree(self):
        for k in (.20,.40,.70):
            with self.subTest(k=k):
                result=m.validation(k)
                self.assertGreater(abs(result["gamma_dB_unit"]["analytic"]),.005)
                self.assertLess(max(result["relative_errors"].values()),5e-4)
                self.assertTrue(result["all_three_nonzero_methods_agree"])

    def test_projector_and_resolvent_local_F_nonzero(self):
        k=.42
        a=m.analytic_F(k,0)
        r=m.F_resolvent(k,0)
        p=m.F_projector(k,0)
        w=m.F_wilson(k,0)
        self.assertGreater(abs(a),.1)
        self.assertAlmostEqual(a,r,places=6)
        self.assertAlmostEqual(a,p,places=6)
        self.assertLess(abs(w-a),5e-6)

    def test_four_independent_wavefunction_phases_do_not_change_wilson_loop(self):
        k=.40
        x=m.wilson_phase(k,0)
        y=m.wilson_phase(k,0,phase_shift=[.7,2.,-.54,1.3])
        self.assertAlmostEqual(x,y,places=12)

    def test_sign_errors_cannot_pass_zero_only_controls(self):
        positive=m.analytic_Gamma(.4)
        self.assertGreater(positive,0.01)
        negative=-m.gamma_centered(m.F_resolvent,.4)
        self.assertLess(negative,0)
        self.assertGreater(abs(negative-positive),0.1)

    def test_explicit_nonmaterial_scientific_boundary(self):
        r=m.validation()
        self.assertEqual(r["new_PBS_jobs"],0)
        self.assertIsNone(r["physical_nmee"])
        self.assertIn("NOT_material".lower(),r["kind"].lower())


if __name__=="__main__":unittest.main()
