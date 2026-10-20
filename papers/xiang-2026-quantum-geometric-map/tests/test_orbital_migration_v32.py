"""Migrated v30 orbital scientific gates. Disagreements must remain visible."""
import unittest
from qgmap.orbital_phe_v4 import orbital_phe_eq56,periodic_data
from qgmap.orbital_mnhe_v5 import orbital_mnhe_2d,orbital_mnhe_3d

class OrbitalGates(unittest.TestCase):
    def test_finite_field_eq56_first_line(self):
        d=orbital_phe_eq56(N=61)
        self.assertLess(abs(d["finite_difference_minus_first"]),5e-7)
        self.assertGreater(abs(d["printed_line_discrepancy"]),1e-3)
        self.assertAlmostEqual(d["printed_line_discrepancy"],-0.03880300933996192,places=10)
    def test_mnhe_2d_sbc(self):
        d=orbital_mnhe_2d(N=61)
        self.assertLess(abs(d["eq52_minus_eq48"]),1e-12)
        self.assertLess(abs(d["finite_field_minus_eq48"]),1e-6)
    def test_mnhe_3d_factor2_and_field(self):
        d=orbital_mnhe_3d(N=17)
        self.assertLess(abs(d["eq49_derived_minus_eq48"]),1e-6)
        self.assertLess(abs(d["finite_field_minus_eq48"]),1e-6)
        self.assertGreater(abs(d["eq52_printed_minus_eq48"]),1e-7)
    def test_sign_and_rejection(self):
        a=orbital_phe_eq56(N=31,a=0,b=1)
        b=orbital_phe_eq56(N=31,a=1,b=0)
        self.assertLess(abs(a["eq56_first_line"]-b["eq56_first_line"]),1e-11)
        with self.assertRaises(ValueError):orbital_phe_eq56(N=19,B_step=0)
        with self.assertRaises(ValueError):orbital_mnhe_3d(N=8)
if __name__=="__main__": unittest.main()
