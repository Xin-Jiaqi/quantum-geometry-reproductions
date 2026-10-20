"""Machine-readable API and scientific-status safety checks for advanced modules."""
import json
import unittest
from qgmap.spin_multiband_v5 import four_band_audit
from qgmap.orbital_phe_v4 import orbital_phe_eq56
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d

class AdvancedContract(unittest.TestCase):
    def test_comparable_same_model_when_repeat(self):
        a=four_band_audit(N=19);b=four_band_audit(N=19)
        self.assertEqual(a["direct_SM29"],b["direct_SM29"])
    def test_all_four_are_json_serializable(self):
        cases=(four_band_audit(N=11),orbital_phe_eq56(N=19),
               orbital_mnhe_2d(N=19),orbital_mnhe_3d(N=9))
        for c in cases:
            self.assertEqual(json.loads(json.dumps(c,allow_nan=False))["protocol"],c["protocol"])
    def test_distinct_physics_contracts(self):
        a=four_band_audit(N=11);b=orbital_phe_eq56(N=19)
        c=orbital_mnhe_2d(N=19);d=orbital_mnhe_3d(N=9)
        self.assertEqual(a["n_band"],4)
        self.assertIn("eq56_first_line",b)
        self.assertIn("eq52_squared_curvature",c)
        self.assertIn("eq52_printed_literal",d)
        self.assertNotEqual(a["units"],b["normalization"])
if __name__=="__main__": unittest.main()
