"""v71 frozen numerical provenance checks independent regenerated observables."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V71Frozen(unittest.TestCase):
    def test_frozen_evidence_independent_regeneration(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        path=root/'scripts/reproduce_v71.py'
        spec=importlib.util.spec_from_file_location('v71_gen',path)
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        expected=json.loads((root/'results/frozen/v71_source_observable_bridge.json').read_text())
        actual=m.generate()
        def check(a,b,p='root'):
            if isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)) and not isinstance(b,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),p)
                self.assertLess(abs(a-b),5e-11,p)
            elif isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),p)
                for k in a:check(a[k],b[k],p+'.'+str(k))
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),p)
                for i,(x,y) in enumerate(zip(a,b)):check(x,y,p+f'[{i}]')
            else:self.assertEqual(a,b,p)
        check(actual,expected)

if __name__=='__main__': unittest.main()
