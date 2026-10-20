"""v40 frozen numerical reference and provenance gate, intentionally partial science."""
import hashlib
import json
import unittest
from pathlib import Path
from qgmap.current_vertex_v40 import current_audit, exact_disorder_ra_torus
ROOT=Path(__file__).resolve().parents[1]

class V40Frozen(unittest.TestCase):
    def test_reference_scba_and_finite_torus(self):
        ref=json.loads((ROOT/'results/frozen/v40_reference.json').read_text())
        self.assertEqual(ref['science_status'],'PARTIAL')
        d=current_audit(mesh=21,energy=1.1,W=.12,eta=.12)
        saved=ref['SCBA']
        for k in ('ra_antisymmetric_dressed','ra_antisymmetric_bare','ra_vertex_effect',
                  'ra_density_ward_residual','four_channel_hall_antisym_abs'):
            with self.subTest(k=k): self.assertAlmostEqual(d[k],saved[k],delta=1e-9)
        exact=exact_disorder_ra_torus(L=4,W=.12,mu=1.1,eta=.12,seeds=6)
        self.assertEqual(exact['samples'],ref['finite_torus']['samples'])
    def test_byte_level_manifest(self):
        manifest=json.loads((ROOT/'results/frozen/v40_manifest.json').read_text())
        self.assertGreaterEqual(len(manifest['files']),5)
        for path,expected in manifest['files'].items():
            data=(ROOT/path).read_bytes()
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256(data).hexdigest(),expected['sha256'])
                self.assertEqual(hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest(),expected['git_blob_sha1'])
if __name__=='__main__':unittest.main()
