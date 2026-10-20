"""Stage-7 source-to-result provenance is verified in GitHub CI."""
from pathlib import Path
import hashlib,json,unittest
ROOT=Path(__file__).resolve().parents[1]

class Stage7Manifest(unittest.TestCase):
    def test_all_frozen_source_and_data_sha256(self):
        m=json.loads((ROOT/'results/STAGE7_MANIFEST.json').read_text())
        self.assertEqual(m['stage'],7)
        self.assertEqual(len(m['sha256']),19)
        for relative,digest in m['sha256'].items():
            with self.subTest(file=relative):
                path=ROOT/relative
                self.assertTrue(path.is_file())
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),digest)

if __name__=='__main__':unittest.main()
