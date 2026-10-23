"""Stage8b: candidate ideal-structure full-motif symmetry regression, no DFT."""
import importlib.util
from pathlib import Path
import unittest

PATH=(Path(__file__).resolve().parents[1]/
      "material-screening"/"ideal_stack_symmetry.py")
spec=importlib.util.spec_from_file_location("stage8b_ideal_stack",PATH)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class IdealStructureSymmetryTests(unittest.TestCase):
    def test_AA_reference_D3h_forbids_gamma_xzx(self):
        for form,row in module.audit().items():
            aa=row["stacks"]["AA"]
            self.assertEqual(aa["order"],12)
            self.assertTrue(aa["C3"])
            self.assertTrue(aa["mz"])
            self.assertFalse(aa["inversion"])
            self.assertFalse(aa["allowed_Gamma_xzx"])

    def test_AB_BA_C3v_allow_and_are_mirror_domains(self):
        for form,row in module.audit().items():
            self.assertTrue(row["AB_to_BA_by_mz"])
            for name in ("AB","BA"):
                state=row["stacks"][name]
                self.assertEqual(state["order"],6)
                self.assertTrue(state["C3"])
                self.assertFalse(state["mz"])
                self.assertFalse(state["inversion"])
                self.assertTrue(state["allowed_Gamma_xzx"])

    def test_generic_slide_has_C1(self):
        for form,row in module.audit().items():
            state=row["stacks"]["generic"]
            self.assertEqual(state["order"],1)
            self.assertTrue(state["allowed_Gamma_xzx"])
            self.assertFalse(state["C3"])


if __name__=="__main__":
    unittest.main()
