"""v76: independent frozen output is regenerated with all real calculations."""
import json,importlib.util,pathlib,unittest
import numpy as np
class V76Frozen(unittest.TestCase):
    def test_regenerate_all_numerical_evidence(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('v76_gen',root/'scripts/reproduce_v76.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        actual=module.generate()
        expected=json.loads((root/'results/frozen/v76_conserving_two_band_collision.json').read_text())
        def compare(a,b,key='root'):
            if isinstance(a,(float,int)) and not isinstance(a,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),key)
                self.assertLess(abs(a-b),2e-8,key)
            elif isinstance(a,dict):
                self.assertEqual(set(a),set(b),key)
                for k in a:compare(a[k],b[k],key+'.'+str(k))
            elif isinstance(a,(list,tuple)):
                self.assertEqual(len(a),len(b),key)
                for i,(x,y) in enumerate(zip(a,b)):compare(x,y,key+str(i))
            else:self.assertEqual(a,b,key)
        compare(actual,expected)
if __name__=='__main__':unittest.main()
