"""Regenerate v67 microscopic coherent lead and Ward frozen evidence in CI."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V67Frozen(unittest.TestCase):
    def test_reservoir_ward_frozen_regeneration(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        modspec=importlib.util.spec_from_file_location('v67_freeze',root/'scripts/reproduce_v67.py')
        module=importlib.util.module_from_spec(modspec);modspec.loader.exec_module(module)
        actual=module.create()
        expected=json.loads((root/'results/frozen/v67_conserving_reservoir_ward.json').read_text())
        def check(x,y,path='root'):
            if type(x) in (int,float) and type(y) in (int,float):
                self.assertTrue(np.isfinite(x) and np.isfinite(y),path)
                # 16-corner mixed derivative has floating subtractive cancellation;
                # analytic Caroli/gauge results retain ~1e-10 tolerance.
                tol=2e-7 if 'fd16_' in path else 1e-10
                self.assertLess(abs(x-y),tol,path)
                return
            self.assertIs(type(x),type(y),path)
            if isinstance(x,dict):
                self.assertEqual(x.keys(),y.keys(),path)
                for k in x: check(x[k],y[k],path+'.'+k)
            elif isinstance(x,list):
                self.assertEqual(len(x),len(y),path)
                for i,(a,b) in enumerate(zip(x,y)):check(a,b,path+str(i))
            else:self.assertEqual(x,y,path)
        check(actual,expected)

if __name__=='__main__':unittest.main()
