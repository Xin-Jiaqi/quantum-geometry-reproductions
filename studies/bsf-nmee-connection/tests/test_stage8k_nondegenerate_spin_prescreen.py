"""Stage8k: zero-job nondegenerate spin evidence checks from original Stage8i JSON."""
import importlib.util
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/"hpc-stage8i"/"stage8k_existing_spin_prescreen.py"
spec=importlib.util.spec_from_file_location("stage8k_existing_spin_prescreen",path)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class NondegenerateSpinEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out=mod.classify()
        cls.cases={(x["job"],x["correct_location"],x["band"]):x
                   for x in cls.out["positive_proxy_candidates"]}

    def test_AB_near_M_non_degenerate_proxy_survives(self):
        x=self.cases[("soc_AB_s0.00","GAMMA_M_NEAR_M__NOT_EXACT_M",35)]
        self.assertAlmostEqual(x["S_perp_pseudized"],.188696,places=4)
        self.assertAlmostEqual(x["min_adjacent_band_separation_meV"],10.29,places=2)
        self.assertFalse(x["verified_physical_spin_matrix"])
        self.assertFalse(x["verified_Sz_nonconservation"])

    def test_AB_Gamma_K_more_isolated_valence_candidate(self):
        x=self.cases[("soc_AB_s0.00","GAMMA_K_INTERIOR__NOT_TRUE_BAND_EDGE",35)]
        self.assertAlmostEqual(x["S_perp_pseudized"],.03360,places=3)
        self.assertAlmostEqual(x["min_adjacent_band_separation_meV"],71.25,places=2)

    def test_original_falsely_named_Kprime_correctly_reinterpreted(self):
        x=self.cases[("soc_AB_s0.00","GAMMA_M_TWO_THIRDS__NOT_K_PRIME",35)]
        self.assertAlmostEqual(x["S_perp_pseudized"],.1678,places=3)
        self.assertAlmostEqual(x["min_adjacent_band_separation_meV"],7.84,places=2)
        self.assertEqual(mod.normalized_label("Kp",[-1/3,-1/3,0]),"TRUE_K_PRIME")

    def test_glide_midpoint_same_false_Kprime_gauge_trap_excluded(self):
        self.assertNotIn(("soc_GLIDE_s0.50","GAMMA_M_TWO_THIRDS__NOT_K_PRIME",35),self.cases)
        matching=[x for x in self.out["matched_k_band_index_comparison"]
                  if x["location"]=="GAMMA_M_TWO_THIRDS__NOT_K_PRIME"
                  and x["band_index"]==35]
        self.assertEqual(len(matching),1)
        self.assertEqual(matching[0]["midpoint_min_adjacent_gap_meV"],0)

    def test_stack_contrast_is_reported_as_orbital_unverified(self):
        comparison=[x for x in self.out["matched_k_band_index_comparison"]
                    if x["location"]=="GAMMA_M_NEAR_M__NOT_EXACT_M"
                    and x["band_index"]==35]
        self.assertEqual(len(comparison),1)
        self.assertGreater(comparison[0]["AB_perp_spin_proxy"],.18)
        self.assertLess(comparison[0]["midpoint_perp_spin_proxy"],1e-4)
        self.assertTrue(comparison[0]["cross_stack_same_band_index_is_not_proven_same_orbital"])

    def test_derived_report_matches_source_and_uses_only_80_sampled_states(self):
        import json
        path=mod.HERE/"STAGE8K_NONDEGENERATE_SPIN_CANDIDATES_20261022.json"
        derived=json.loads(path.read_text())
        self.assertEqual(derived["scan"]["archived_samples_processed"],80)
        self.assertEqual(self.out["sampled_states_count"],80)
        self.assertEqual(derived["scan"]["total_ab_nondegenerate_proxy_candidates"],
                         sum(x["job"]=="soc_AB_s0.00" for x in
                             self.out["positive_proxy_candidates"]))
        self.assertEqual(derived["scan"]["total_glide_nondegenerate_proxy_candidates"],
                         sum(x["job"]=="soc_GLIDE_s0.50" for x in
                             self.out["positive_proxy_candidates"]))
        self.assertEqual(derived["scan"]["total_glide_nondegenerate_proxy_candidates"],0)
        self.assertTrue(any(x["location"]=="GAMMA_K_INTERIOR__NOT_TRUE_BAND_EDGE"
                            for x in derived["priority_examples"]))

    def test_no_new_hpc_jobs_and_no_claim_of_material_gamma(self):
        self.assertEqual(self.out["max_additional_PBS_jobs"],0)
        self.assertEqual(self.out["static_material_NMEE"],"NOT_COMPUTED")
        self.assertIn("NOT_VERIFIED",self.out["spin_mixing_gate"])
        self.assertTrue(all(not x["verified_physical_spin_matrix"]
                            for x in self.out["positive_proxy_candidates"]))

if __name__=="__main__":unittest.main()
