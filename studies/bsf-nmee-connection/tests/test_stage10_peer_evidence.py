"""Stage10 peer-review evidence bookkeeping checks; scientific claims remain cited in report."""
import json
from pathlib import Path
import unittest

DOC=Path(__file__).resolve().parents[1]/"docs"
class Stage10PeerEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((DOC/"STAGE10_PEER_MATRIX.json").read_text())
        cls.text=(DOC/"STAGE10_PEER_LANDSCAPE_20261021.md").read_text()

    def test_sources_have_unique_identifiers_and_traceable_primary_urls(self):
        sources=self.data["matrix"]
        self.assertGreaterEqual(len(sources),14)
        ids=[x["id"] for x in sources]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertTrue(all(x["evidence"] and x["title"] and x["doi_or_url"] for x in sources))
        for x in sources:
            if x["doi_or_url"].startswith("10."):
                self.assertIn("https://doi.org/"+x["doi_or_url"],self.text)

    def test_critical_adjacent_2026_competitors_are_not_missing(self):
        byid={x["id"]:x for x in self.data["matrix"]}
        for key in ("P01","P02","P03","P04","P05","P07","P09"):
            self.assertIn(key,byid)
            self.assertIn(byid[key]["title"],self.text)
        self.assertIn("10.1103/2gyh-y52n",byid["P03"]["doi_or_url"])
        self.assertIn("B_arxiv_preprint",byid["P09"]["evidence"])

    def test_does_not_claim_unproven_material_nmee_or_static_dc_equality(self):
        status=self.data["publication_hypothesis_status"]
        self.assertEqual(status["first_intrinsic_insulator_NMEE"],"PREEMPTED")
        self.assertEqual(status["first_layer_hidden_response"],"PREEMPTED")
        self.assertIn("UNCALCULATED",status["remote_band_spin_mixing_detailed_mechanism"])
        self.assertEqual(self.data["full_formula_sm_review_status"],"NOT_COMPLETED")
        self.assertIn("不是经证实等于",self.text)
        self.assertIn("尚无实际非零数值",self.text)

if __name__=="__main__":
    unittest.main()
