"""Stage8l: guard against off-by-one band indexing in published audit JSON."""
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/"hpc-stage8i"/"stage8l_crosscheck_band_ledger.py"
spec=importlib.util.spec_from_file_location("stage8l_crosscheck_band_ledger",PATH)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class Stage8lBandEvidenceConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=mod.verify()
    def test_new_band35_label_actually_matches_old_band36(self):
        self.assertEqual(len(self.report["unexpected_band35_label_issues"]),7)
        self.assertEqual(len(self.report["all_records"]),8)
        self.assertTrue(all(x["independent_saved_S_grouped_corresponds_to_band"]==36
                            for x in self.report["unexpected_band35_label_issues"]))
        ab=[x for x in self.report["unexpected_band35_label_issues"]
             if x["job"]=="soc_AB_s0.00"]
        self.assertEqual(len(ab),4)
        for x in ab:
            self.assertTrue(x["saved_S_perp"]>=0)
    def test_degenerate_glide_pair_is_not_arbitrarily_relabeled(self):
        amb=[x for x in self.report["all_records"]
             if x["degenerate_energy_only_cannot_disambiguate"]]
        self.assertEqual(len(amb),1)
        self.assertEqual(amb[0]["job"],"soc_GLIDE_s0.50")
        self.assertIn("OLD_FALSE_Kp",amb[0]["location"])
    def test_no_vasp_and_no_static_nmee_claims(self):
        self.assertEqual(self.report["new_vasp_jobs"],0)
        self.assertFalse(self.report["raw_wavecar_examined_by_this_check"])
        self.assertEqual(self.report["Kubo_static_NMEE"],"NOT_COMPUTED")
        self.assertEqual(self.report["PAW_physical_spin"],"NOT_QUANTITATIVELY_VALIDATED")

if __name__=="__main__":unittest.main()
