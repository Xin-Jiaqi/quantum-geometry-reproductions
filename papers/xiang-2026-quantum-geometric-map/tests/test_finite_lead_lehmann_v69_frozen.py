"""Independent fast regeneration of frozen v69 finite-lead Lehmann evidence."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V69Frozen(unittest.TestCase):
    def test_quick_lehmann_frozen_evidence_regenerates(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('v69_reproduce',root/'scripts/reproduce_v69.py')
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        expected=json.loads((root/'results/frozen/v69_finite_lead_lehmann.json').read_text())
        actual=mod.create(False)
        def rec(a,b,where='root'):
            if isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)) and not isinstance(b,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),where)
                self.assertLess(abs(a-b),2e-9,where);return
            self.assertIs(type(a),type(b),where)
            if isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),where)
                for k in a:rec(a[k],b[k],where+'.'+k)
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),where)
                for i,(x,y) in enumerate(zip(a,b)):rec(x,y,where+f'[{i}]')
            else:self.assertEqual(a,b,where)
        rec(actual['quick'],expected['quick'])
        rec(actual['keldysh_reference'],expected['keldysh_reference'])

if __name__=='__main__':unittest.main()
