import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from integration.stage1_handoff import inspect_package,REQUIRED

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

class HandoffTest(unittest.TestCase):
 def fixture(self,root):
    w=root/'L1';w.mkdir()
    (w/'POSCAR').write_text('Synthetic\n1\n1 0 0\n0 1 0\n0 0 20\nH\n1\nDirect\n0 0 0\n')
    (w/'wannier90_hr.dat').write_text('synthetic\n2\n1\n1\n0 0 0 1 1 -1 0\n0 0 0 2 1 0 0\n0 0 0 1 2 0 0\n0 0 0 2 2 1 0\n')
    (w/'wannier90-wcenter.dat').write_text('X 0 0 0\nX 1 0 0\n')
    (w/'quality.json').write_text(json.dumps({'passed':True,'num_wann':2,'nocc_model':1,'soc':True}))
    m={'num_wann':2,'nocc_model':1,'physics':{'soc':True,'dimensionality':2},
       'structures':{'L1':{'files':{n:sha(w/n) for n in REQUIRED}}}}
    (root/'package-manifest.json').write_text(json.dumps(m))

 def test_synthetic_handoff_preflight_fail_closed(self):
    with tempfile.TemporaryDirectory() as d:
       root=Path(d);self.fixture(root)
       result=inspect_package(root)
       self.assertEqual(result['status'],'BLOCKED_PHYSICAL_OPERATORS')
       self.assertEqual(result['num_wann'],2)
       self.assertTrue(result['sha256_verified'])

 def test_bad_hash_rejected(self):
    with tempfile.TemporaryDirectory() as d:
       root=Path(d);self.fixture(root)
       (root/'L1'/'quality.json').write_text('{}')
       with self.assertRaisesRegex(ValueError,'SHA256'):inspect_package(root)

 def test_inconsistent_soc_rejected(self):
    with tempfile.TemporaryDirectory() as d:
       root=Path(d);self.fixture(root)
       f=root/'L1'/'quality.json';q=json.loads(f.read_text());q['soc']=False
       f.write_text(json.dumps(q))
       path=root/'package-manifest.json';m=json.loads(path.read_text());m['structures']['L1']['files']['quality.json']=sha(f)
       path.write_text(json.dumps(m))
       with self.assertRaisesRegex(ValueError,'SOC conflict'):inspect_package(root)

if __name__=='__main__':unittest.main()
