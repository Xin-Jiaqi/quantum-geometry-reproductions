"""Recompute v74 full frozen pole-regulator data via generator and compare."""
import unittest,importlib.util,pathlib,json
import numpy as np

class V74Frozen(unittest.TestCase):
    def test_independent_reproduction_matches_numeric_frozen(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        modspec=importlib.util.spec_from_file_location('v74_gen',root/'scripts/reproduce_v74.py')
        mod=importlib.util.module_from_spec(modspec);modspec.loader.exec_module(mod)
        actual=mod.generate()
        expected=json.loads((root/'results/frozen/v74_selected_electric_pole.json').read_text())
        def rec(a,b,p='root'):
            if isinstance(a,(int,float)) and not isinstance(a,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),p)
                self.assertLess(abs(a-b),1e-10,p)
            elif isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),p)
                for k in a:rec(a[k],b[k],p+'.'+k)
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),p)
                for i,(x,y) in enumerate(zip(a,b)):rec(x,y,p+f'[{i}]')
            else:self.assertEqual(a,b,p)
        rec(actual,expected)
