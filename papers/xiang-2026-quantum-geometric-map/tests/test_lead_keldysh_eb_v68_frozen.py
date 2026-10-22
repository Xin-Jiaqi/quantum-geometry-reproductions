"""Frozen independent recomputation v68, tolerance separated by quadrature type."""
import importlib.util,json,pathlib,unittest
import numpy as np

class V68Frozen(unittest.TestCase):
    def test_full_noninteracting_keldysh_frozen_evidence(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        sp=importlib.util.spec_from_file_location('regenerate_v68',root/'scripts/reproduce_v68.py')
        mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
        got=mod.create()
        expected=json.loads((root/'results/frozen/v68_same_device_keldysh_eb.json').read_text())
        def rec(x,y,where='root'):
            if isinstance(x,(int,float)) and not isinstance(x,bool) and isinstance(y,(int,float)) and not isinstance(y,bool):
                self.assertTrue(np.isfinite(x) and np.isfinite(y),where)
                # Five-point or finite gate quadrature has independent floating-error
                # scale; do not weaken analytic or microscopic Ward identities.
                tol=2e-7 if any(k in where for k in ('finite_gate','five_point','convergence.384','convergence.64')) else 5e-9
                self.assertLess(abs(x-y),tol,where)
                return
            self.assertIs(type(x),type(y),where)
            if isinstance(x,dict):
                self.assertEqual(x.keys(),y.keys(),where)
                for k in x:rec(x[k],y[k],where+'.'+k)
            elif isinstance(x,list):
                self.assertEqual(len(x),len(y),where)
                for i,(a,b) in enumerate(zip(x,y)):rec(a,b,where+str(i))
            else:self.assertEqual(x,y,where)
        rec(got,expected)

if __name__=='__main__':unittest.main()
