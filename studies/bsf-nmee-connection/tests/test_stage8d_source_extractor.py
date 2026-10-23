"""Source-only BiDB screening gate, requiring no ASE installation."""
import unittest
from pathlib import Path
import importlib.util,sqlite3
from tempfile import TemporaryDirectory

MODULE=Path(__file__).resolve().parents[1]/"material-screening"/"extract_bidb_targets.py"
spec=importlib.util.spec_from_file_location("stage8d_extract",MODULE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

class BiDBExtractionGates(unittest.TestCase):
    def test_exact_six_atom_target_only(self):
        self.assertEqual(mod.identify(["W"]*2+["Se"]*4),"WSe2")
        self.assertEqual(mod.identify(["Mo"]*2+["S"]*4),"MoS2")
        self.assertIsNone(mod.identify(["W"]+["Se"]*2))
        self.assertIsNone(mod.identify(["W"]*3+["Se"]*6))

    def test_unknown_sqlite_not_promoted_to_structure(self):
        with TemporaryDirectory() as folder:
            path=Path(folder)/"empty.db"
            sqlite3.connect(path).close()
            self.assertFalse(mod.ase_sqlite_shape(path))
            self.assertEqual(len(mod.sha256_file(path)),64)

    def test_schema_is_gate_not_relaxation_claim(self):
        with TemporaryDirectory() as folder:
            path=Path(folder)/"synthetic.db"
            con=sqlite3.connect(path)
            con.execute("CREATE TABLE systems (id INTEGER)")
            con.execute("CREATE TABLE species (id INTEGER)")
            con.commit();con.close()
            self.assertTrue(mod.ase_sqlite_shape(path))

if __name__=="__main__":
    unittest.main()
