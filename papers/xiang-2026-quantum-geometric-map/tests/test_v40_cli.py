"""Standalone v40 command JSON must remain explicitly scientifically PARTIAL."""
import json
from pathlib import Path
import subprocess
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]

class CLI(unittest.TestCase):
    def test_local_cli(self):
        p=subprocess.run([sys.executable,str(ROOT/'scripts/reproduce_current_vertex_v40.py'),
                          '--mesh','17','--L','3','--seeds','2'],cwd=ROOT,
                          text=True,capture_output=True,timeout=40)
        self.assertEqual(p.returncode,0,p.stderr)
        d=json.loads(p.stdout)
        self.assertEqual(d['science_status'],'PARTIAL')
        self.assertIn('NOT a physical Hall conductivity',d['SCBA']['warnings'][0])
        self.assertLess(d['SCBA']['four_channel_hall_antisym_abs'],1e-10)
if __name__=='__main__':unittest.main()
