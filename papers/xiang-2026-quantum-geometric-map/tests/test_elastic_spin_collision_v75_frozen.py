"""Independent regenerate source-sign and elastic collision numerical frozen evidence."""
import json,pathlib,importlib.util,unittest
import numpy as np

class V75Frozen(unittest.TestCase):
 def test_frozen_regeneration(self):
  root=pathlib.Path(__file__).resolve().parents[1]
  spec=importlib.util.spec_from_file_location('v75_reproduce',root/'scripts/reproduce_v75.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  actual=module.generate();expected=json.loads((root/'results/frozen/v75_elastic_collision_and_source_sign.json').read_text())
  def rec(a,b,ref='root'):
   if isinstance(a,(int,float)) and not isinstance(a,bool):
    self.assertTrue(np.isfinite(a) and np.isfinite(b),ref)
    self.assertLess(abs(a-b),1e-10,ref)
   elif isinstance(a,dict):
    self.assertEqual(a.keys(),b.keys(),ref)
    for k in a:rec(a[k],b[k],ref+'.'+k)
   elif isinstance(a,list):
    self.assertEqual(len(a),len(b),ref)
    for i,(x,y) in enumerate(zip(a,b)):rec(x,y,ref+str(i))
   else:self.assertEqual(a,b,ref)
  rec(actual,expected)
if __name__=='__main__':unittest.main()
