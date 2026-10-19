import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bsf_nmee.connection import _validate_operation, transform_layer_responses, Connection, response_summary


class TestCoordinateAndTensorContract(unittest.TestCase):
    def setUp(self):
        r=np.diag([1.,1.,-1.])
        self.c=Connection((1/3,2/3),(2/3,1/3),'m_001',True,r,r,
                          {'top':'bottom','bottom':'top'},'P')
    def test_mirror_z_total_odd_difference_even(self):
        a=np.zeros((3,3,3));b=np.zeros_like(a)
        a[0,2,0]=2.;b[0,2,0]=.5
        mapped=transform_layer_responses({'top':a,'bottom':b},self.c)
        old=response_summary({'top':a,'bottom':b})
        new=response_summary(mapped)
        self.assertAlmostEqual(new['total'][0,2,0],-old['total'][0,2,0])
        self.assertAlmostEqual(new['layer_difference'][0,2,0],old['layer_difference'][0,2,0])
        self.assertAlmostEqual(mapped['top'][0,2,0],-.5)
        self.assertAlmostEqual(mapped['bottom'][0,2,0],-2.)
    def test_mirror_involution(self):
        rng=np.random.default_rng(10)
        a={'top':rng.normal(size=(3,3,3)),'bottom':rng.normal(size=(3,3,3))}
        twice=transform_layer_responses(transform_layer_responses(a,self.c),self.c)
        for k in a:self.assertLess(np.max(abs(twice[k]-a[k])),1e-13)
    def test_identity_and_inversion_transform(self):
        for R in (np.eye(3),-np.eye(3)):
            a=np.arange(27.).reshape(3,3,3)
            c=Connection((.2,.3),(.2,.3),'operation',R[2,2]<0,R,R,
                        {'top':'bottom','bottom':'top'} if R[2,2]<0 else {'top':'top','bottom':'bottom'},'P')
            out=transform_layer_responses({'top':a,'bottom':-a},c)
            expected=R[0,0]*R[2,2]*R[0,0]*a[0,2,0]
            dest=c.layer_permutation['top']
            self.assertAlmostEqual(out[dest][0,2,0],expected)
    def test_reject_unsafe_rotation_and_missing_basis(self):
        with self.assertRaises(TypeError):_validate_operation(np.diag([1.,1.,-1.]))
        with self.assertRaises(ValueError):_validate_operation(SimpleNamespace(name='bad',matrix_fractional=np.diag([1,1,-1]),matrix_cartesian=np.array([[1,0,.1],[0,1,0],[0,0,-1]])))
    def test_invalid_layer_or_tensor(self):
        with self.assertRaises(ValueError):transform_layer_responses({'top':np.zeros((3,3,3))},self.c)
        with self.assertRaises(ValueError):transform_layer_responses({'top':np.zeros((2,2,2)), 'bottom':np.zeros((3,3,3))},self.c)

if __name__=='__main__':unittest.main()
