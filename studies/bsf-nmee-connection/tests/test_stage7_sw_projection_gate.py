"""Narrow regression of fixed-projector SW applicability and symmetry."""
import unittest
import numpy as np
from bsf_nmee.sw_projection_diagnostics import sw_gate


class Stage7SWGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = sw_gate(u_values=(.1,), tunneling_values=(.005, .02),
                           q_x_values=(0., .1))

    def test_sample_gap_and_transformations(self):
        self.assertEqual(len(self.rows), 8)
        self.assertTrue(all(np.isfinite(r["sw_diagnostic_ratio"]) for r in self.rows))
        self.assertTrue(all(r["mu_distance_meV"] > 0 for r in self.rows))
        self.assertLess(max(r["mirror_error_eV"] for r in self.rows), 1e-10)
        self.assertLess(max(r["TR_error_eV"] for r in self.rows), 1e-10)

    def test_fourband_fixed_projection_not_perturbative_at_20mev(self):
        v = next(x for x in self.rows if x["valley"] == 0
                 and x["t_perp_meV"] == 20 and x["q_x_inv_nm"] == 0)
        self.assertGreater(v["sw_diagnostic_ratio"], 1)
        self.assertLess(v["smallest_principal_cos_squared"], 0.01)

    def test_valleys_have_same_spectra_and_diagnostics(self):
        for r in self.rows:
            match = next(x for x in self.rows if x["valley"] == 1-r["valley"]
                         and x["u"] == r["u"]
                         and x["q_x_inv_nm"] == r["q_x_inv_nm"]
                         and x["t_perp_meV"] == r["t_perp_meV"])
            for field in ("full_middle_gap_meV", "mu_distance_meV",
                          "sw_diagnostic_ratio"):
                self.assertAlmostEqual(r[field], match[field], places=7)


if __name__ == "__main__":
    unittest.main()
