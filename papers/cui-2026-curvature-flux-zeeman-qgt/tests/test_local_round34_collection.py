"""Round34 local harness receipt is offline, bounded and deliberately unauthenticated."""
import tempfile,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from collect_round34_local_sources import collect

class LocalHarnessReadOnlyReceipt(unittest.TestCase):
 def test_complete_fixture_reports_candidate_not_authentic(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);p=root/"seed"
   for suff in (".spn.fmt","_u.mat","_u_dis.mat","_v.mat",
                ".eig","_hr.dat","_r.dat",".win",".wout"):
    (root/(p.name+suff)).write_text("FAKE ASCII CONTENT")
   r=collect([root])
   self.assertEqual(r["found_seed_groups"],1)
   self.assertFalse(r["real_material_KME_verified"])
   self.assertTrue(r["groups"][0]["plausible_ascii_chain"])
   self.assertFalse(r["groups"][0]["physical_spin_source_authenticated"])
   self.assertFalse(r["groups"][0]["can_claim_real_kme"])
   self.assertFalse(r["original_files_uploaded"])

 def test_missing_spin_and_computation_never_fake_provenance(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   (p/"x_hr.dat").write_text("H")
   (p/"x_u.mat").write_text("U")
   r=collect([p])
   self.assertEqual(r["found_seed_groups"],1)
   self.assertFalse(r["groups"][0]["plausible_ascii_chain"])
   self.assertIn("spn_fmt",r["groups"][0]["missing_or_unverified"])
   self.assertFalse(r["computations_submitted"])

 def test_original_ascii_seed_spn_detected_as_formatted(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   (p/"seed.spn").write_text("Created by unit fixture\n4 3\n0.0 0.0\n")
   r=collect([p])
   self.assertIn("spn_fmt",r["groups"][0]["present"])
   self.assertNotIn("spn_binary",r["groups"][0]["present"])

 def test_fortran_binary_spn_is_not_silently_parsed(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   (p/"seed.spn").write_bytes(b"\x00\x80\x90\x00")
   r=collect([p])
   self.assertIn("spn_binary",r["groups"][0]["present"])
   self.assertIn("FORMATTED_SPN_CONVERSION_REQUIRED_OFFLINE",
                 r["groups"][0]["missing_or_unverified"])

 def test_sha_records_not_source_authenticity(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/"fake_u.mat").write_text("obviously fake")
   r=collect([p])
   self.assertEqual(len(r["groups"][0]["present"]["u"]["sha256"]),64)
   self.assertFalse(r["groups"][0]["gauge_and_material_provenance_authenticated"])

 def test_directory_symlink_not_followed(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   (p/"inner").mkdir()
   (p/"inner"/"s.spn.fmt").write_text("inside")
   (p/"shortcut").symlink_to(p/"inner",target_is_directory=True)
   r=collect([p])
   self.assertEqual(r["found_seed_groups"],1)

 def test_home_or_root_directory_safety_rejected(self):
  with self.assertRaisesRegex(ValueError,"not / or entire home"):
   collect(["/"])
  with self.assertRaisesRegex(ValueError,"not / or entire home"):
   collect([Path.home()])

 def test_bounded_scan_fails_nonfatally_without_promoting_material(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)
   for i in range(5):(p/(f"{i}.spn.fmt")).write_text("fake")
   r=collect([p],max_files=2)
   self.assertTrue(r["search_cutoff_reached"])
   self.assertFalse(r["real_material_KME_verified"])

if __name__=="__main__":unittest.main()
