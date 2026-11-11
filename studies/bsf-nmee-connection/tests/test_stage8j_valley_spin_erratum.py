"""Stage8j evidence review: prevent recurrence of impossible time-reversal valley claims."""
import importlib.util
from pathlib import Path
import unittest

SRC=Path(__file__).resolve().parents[1]/"hpc-stage8i"/"stage8j_evidence_check.py"
spec=importlib.util.spec_from_file_location("stage8j_k_label_correction",SRC)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class Stage8jScientificErratumTests(unittest.TestCase):
    def test_mislabeled_kprime_is_interior_gamma_M_line(self):
        r=mod.evaluate()
        self.assertEqual(r["scientific_status"],
            "SOC_SAMPLED_GAP_PASS__SPIN_MIXING_NOT_YET_ESTABLISHED__OUTPUT_COMPLETION_WARNING")
        k=r["genuine_K"]
        true_kp=r["genuine_time_reversed_Kprime"]
        self.assertAlmostEqual(k[0]+true_kp[0],0,places=8)
        self.assertAlmostEqual(k[1]+true_kp[1],0,places=8)
        self.assertLess(r["false_Kp_norm_inverse_A"],r["K_norm_inverse_A"])
        self.assertAlmostEqual(r["false_Kp_is_at_fraction_of_Gamma_M"],2/3)
        self.assertFalse(r["eigenvalue_TR_equality_rechecked_directly_against_original_EIGENVAL"])

    def test_k_band35_36_both_same_spin_not_soc_partner_splitting(self):
        r=mod.evaluate()
        for item in r["jobs"].values():
            self.assertGreater(item["band35_Sz"]*item["band36_Sz"],0)
            self.assertFalse(item["these_bands_form_verified_opposite_spin_pair"])
            self.assertFalse(item["this_is_true_valley_split"])
            self.assertEqual(item["original_report_spin_split_label"],"NOT_SUPPORTED")
        self.assertAlmostEqual(
            r["jobs"]["soc_AB_s0.00"]["band35_band36_energy_separation_meV"],
            87.54,places=2)
        self.assertAlmostEqual(
            r["jobs"]["soc_GLIDE_s0.50"]["band35_band36_energy_separation_meV"],
            35.31,places=2)

    def test_degenerate_spin_rotation_counterexample_preserves_Sz(self):
        r=mod.evaluate()
        c=r["spin_conservation_counterexample"]
        self.assertTrue(c["Sz_commutator_exactly_zero"])
        self.assertEqual(c["degenerate_state_Sx"],0.5)

    def test_no_new_jobs_no_claim_real_gamma(self):
        r=mod.evaluate()
        self.assertEqual(r["new_PBS_jobs"],0)
        self.assertFalse(r["actual_material_Gamma_computed"])

if __name__=="__main__":unittest.main()
