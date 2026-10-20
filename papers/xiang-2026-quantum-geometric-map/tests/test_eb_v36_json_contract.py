"""v36 CLI contract, original v21-v23 model numerical references, protocol labels."""
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FROZEN=json.loads((ROOT/'results/frozen/v36_reference.json').read_text())

def call(name):
    r=subprocess.run([sys.executable,str(ROOT/'scripts/reproduce_advanced.py'),
        '--case',name],capture_output=True,text=True,cwd=ROOT,timeout=90)
    if r.returncode:raise AssertionError(r.stdout+'\n'+r.stderr)
    return json.loads(r.stdout)

def val(x): return complex(x['real'],x['imag'])

class V36Cli(unittest.TestCase):
    def test_v21_analytic_ward_and_gauge(self):
        d=call('eb-finite-q-v21')
        self.assertEqual(d['science_status'],'PARTIAL')
        r=d['result'];f=FROZEN['v21']
        self.assertLess(abs(val(r['liouville'])-val(f['liouville_slope'])),2e-11)
        self.assertLess(r['finite_difference_error'],f['finite_B_difference_tol'])
        self.assertLess(r['mixed_ward_density_abs'],f['mixed_ward_tol'])
        self.assertLess(r['pure_gauge_magnetic_total_abs'],f['gauge_total_tol'])
    def test_v22_ordering_independent(self):
        d=call('eb-two-color-v22');r=d['result'];f=FROZEN['v22']
        for k in ('B_then_E','E_then_B','total'):
            self.assertLess(abs(val(r[k])-val(f[k])),2e-11)
        self.assertLess(r['independent_sylvester_error'],f['independent_error_tol'])
    def test_v23_homogeneous_eq20(self):
        d=call('eb-uniform-v23');r=d['result'];f=FROZEN['v23']
        self.assertLess(abs(val(r['total'])-val(f['total'])),2e-11)
        self.assertLess(abs(val(r['sm20_pole'])-f['sm20_pole']),1e-12)
        self.assertLess(r['independent_algebraic_error'],f['algebraic_error_tol'])
    def test_constraints_and_status(self):
        for name in ('eb-finite-q-v21','eb-two-color-v22','eb-uniform-v23'):
            d=call(name)
            self.assertEqual(d['scientific_status'] if 'scientific_status' in d else d['science_status'],'PARTIAL')
            self.assertIn('limitations',d)
            self.assertIn('result',d)

if __name__=='__main__':unittest.main()