"""v73 independently regenerate every value in frozen two-color evidence."""
import importlib.util,json,pathlib,unittest
import numpy as np
class V73Frozen(unittest.TestCase):
    def test_full_frozen_independent_regeneration(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('v73repro',root/'scripts/reproduce_v73.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        actual=module.generate()
        expected=json.loads((root/'results/frozen/v73_two_color_finite_q_spin.json').read_text())
        def check(a,b,name='root'):
            if isinstance(a,(int,float)) and not isinstance(a,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),name)
                self.assertLess(abs(a-b),2e-10,name)
            elif isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),name)
                for k in a:check(a[k],b[k],name+'.'+str(k))
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),name)
                for i,(x,y) in enumerate(zip(a,b)):check(x,y,name+f'[{i}]')
            else:self.assertEqual(a,b,name)
        check(actual,expected)
if __name__=='__main__':unittest.main()
