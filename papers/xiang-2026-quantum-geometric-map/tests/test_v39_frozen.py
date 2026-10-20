"""Frozen v39 model numbers; intentionally not a complete scientific-acceptance gate."""
import json,subprocess,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class V39Frozen(unittest.TestCase):
    def test_frozen_numbers_and_protocol_label(self):
        data=json.loads((ROOT/'results/frozen/v39_benchmark.json').read_text())
        self.assertEqual(data['scientific_status'],'PARTIAL')
        self.assertAlmostEqual(data['derivative_W']['gamma_minus_identity'],.3295819298185,delta=1e-5)
        self.assertAlmostEqual(data['clean_insulator']['derivative_correction'],-.0266699537731,delta=4e-5)
        self.assertAlmostEqual(data['metallic_control']['derivative_correction'],-.123374195115,delta=5e-5)
        self.assertIn('transverse current ladder',data['metallic_control']['explicit_missing'])
    def test_cli_schema_and_real_calculation(self):
        run=subprocess.run([sys.executable,str(ROOT/'scripts/reproduce_scba_streda_v39.py'),
            '--mesh','17','--quadrature','160'],capture_output=True,text=True,cwd=ROOT,timeout=60)
        self.assertEqual(run.returncode,0,msg=run.stderr)
        d=json.loads(run.stdout)
        self.assertEqual(d['schema'],'qgmap-scba-sea-v39-1')
        self.assertEqual(d['scientific_status'],'PARTIAL')
        self.assertGreater(abs(d['metallic_control']['derivative_correction']),.1)
if __name__=='__main__':unittest.main()
