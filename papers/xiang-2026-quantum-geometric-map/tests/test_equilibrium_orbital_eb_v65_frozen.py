"""v65 frozen-provenance acceptance: independent regeneration, no fitted constants."""
import importlib.util,json,pathlib,unittest
import numpy as np


class V65FrozenAcceptance(unittest.TestCase):
    def test_frozen_independent_reconstruction(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        path=root/'scripts/reproduce_v65.py'
        spec=importlib.util.spec_from_file_location('v65_freeze_reproduce',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        actual=mod.generate()
        expected=json.loads((root/'results/frozen/v65_static_equilibrium_orbital_eb.json').read_text())
        def compare(a,b,label='root'):
            self.assertIs(type(a),type(b),label)
            if isinstance(a,dict):
                self.assertEqual(a.keys(),b.keys(),label)
                for k in a:compare(a[k],b[k],label+'.'+k)
            elif isinstance(a,list):
                self.assertEqual(len(a),len(b),label)
                for k,(u,v) in enumerate(zip(a,b)):compare(u,v,label+'['+str(k)+']')
            elif isinstance(a,float):
                self.assertTrue(np.isfinite(a) and np.isfinite(b),label)
                self.assertLess(abs(a-b),2e-8,label)
            else:self.assertEqual(a,b,label)
        compare(actual,expected)


if __name__=='__main__':unittest.main()
