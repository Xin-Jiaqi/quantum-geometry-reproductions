"""Independently reconstruct all frozen v72 finite-q source and Ward evidence."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V72FrozenAcceptance(unittest.TestCase):
    def test_numerical_provenance_from_independent_rerun(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        spec=importlib.util.spec_from_file_location('v72_rerun',root/'scripts/reproduce_v72.py')
        obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
        expected=json.loads((root/'results/frozen/v72_bulk_finite_q_spin.json').read_text())
        actual=obj.generate()
        def cmp(a,b,label='root'):
            if isinstance(a,(int,float)) and not isinstance(a,bool):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),label)
                self.assertLess(abs(a-b),5e-10,label)
            elif isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),label)
                for k in a:cmp(a[k],b[k],label+'.'+str(k))
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),label)
                for i,(x,y) in enumerate(zip(a,b)):cmp(x,y,label+str(i))
            else:self.assertEqual(a,b,label)
        cmp(actual,expected)

if __name__=='__main__':unittest.main()
