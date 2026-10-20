"""Release checksum gate only for files newly added in stage 5.

No prior frozen results/source file is modified.
"""
from pathlib import Path
import hashlib,json,unittest
ROOT=Path(__file__).resolve().parents[1]

class TestStage5Release(unittest.TestCase):
    def test_checksums(self):
        path=ROOT/'results/STAGE5_MANIFEST.json'
        m=json.loads(path.read_text(encoding='utf8'))
        self.assertEqual(m['schema'],'bsf-nmee-stage5-v1')
        for relative,digest in m['sha256'].items():
            f=ROOT/relative
            self.assertTrue(f.exists(),relative)
            # Normalize repository text EOLs if transferred via GitHub content API.
            normalized=f.read_text(encoding='utf8').replace('\r\n','\n').encode('utf8')
            self.assertEqual(hashlib.sha256(normalized).hexdigest(),digest,relative)

if __name__=='__main__': unittest.main()
