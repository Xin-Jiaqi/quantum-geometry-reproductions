"""v40 advanced JSON entry is executable and guarded by PARTIAL science status."""
import json
import subprocess
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class AdvancedCurrentV40(unittest.TestCase):
    def test_scoped_json(self):
        cmd=[sys.executable,str(ROOT/"scripts/reproduce_advanced.py"),
             "--case","scba-current-vertex-v40","--mesh","17"]
        r=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=60)
        self.assertEqual(r.returncode,0,msg=r.stderr)
        d=json.loads(r.stdout)
        self.assertEqual(d['science_status'],'PARTIAL')
        s=d['result']
        self.assertEqual(s['scientific_status'],'PARTIAL')
        self.assertGreater(abs(s['SCBA']['ra_vertex_effect']),.1)
        self.assertLess(s['SCBA']['four_channel_hall_antisym_abs'],1e-10)
        self.assertIn('NOT full physical Hall',s['physics_scope'])
if __name__=='__main__':unittest.main()
