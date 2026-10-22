"""Independent remote CI regeneration of v66 frozen model evidence."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V66Frozen(unittest.TestCase):
    def test_full_spectral_contact_evidence_regenerates(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('frozen_v66',root/'scripts/reproduce_v66.py')
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        expected=json.loads((root/'results/frozen/v66_dynamic_orbital_mixed_EB.json').read_text())
        actual=mod.create()
        def rec(a,b,where='root'):
            if isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)) and not isinstance(b,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),where)
                self.assertLess(abs(a-b),1e-10,where)
                return
            self.assertIs(type(a),type(b),where)
            if isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),where)
                for key in a:rec(a[key],b[key],where+'.'+key)
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),where)
                for i,(x,y) in enumerate(zip(a,b)):rec(x,y,where+str(i))
            elif isinstance(a,float):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),where)
                self.assertLess(abs(a-b),1e-10,where)
            else:self.assertEqual(a,b,where)
        rec(actual,expected)

if __name__=='__main__':unittest.main()
