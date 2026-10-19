"""These tests require the ACTUAL public group-theory toolkit (no replacement)."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bsf_nmee import connect_interface, transform_layer_responses

try:
    from group_theory_operations import point_group_operations
except ImportError:
    point_group_operations=None

@unittest.skipUnless(point_group_operations is not None, 'Install group-theory-operations-toolkit v0.18.0')
class ActualToolkitTests(unittest.TestCase):
    def setUp(self):
        self.ops={op.name:op for op in point_group_operations('6/mmm')}
    def test_bn_partner_mirror_z(self):
        c=connect_interface((1/3,2/3),self.ops['m_001'],expected_target=(2/3,1/3))
        self.assertTrue(c.layer_exchanging)
        self.assertEqual(c.layer_permutation,{'top':'bottom','bottom':'top'})
        self.assertLess(abs(c.rotation_cartesian[2,2]+1),1e-14)
    def test_partner_mirror_involution(self):
        c1=connect_interface((1/3,2/3),self.ops['m_001'])
        c2=connect_interface(c1.target_shift,self.ops['m_001'],expected_target=(1/3,2/3))
        self.assertTrue(np.allclose(c2.target_shift,(1/3,2/3)))
        self.assertTrue(c2.layer_exchanging)
    def test_identity_and_expected_target_rejection(self):
        ident=connect_interface((1/3,2/3),self.ops['1'],expected_target=(1/3,2/3))
        self.assertFalse(ident.layer_exchanging)
        with self.assertRaises(ValueError):connect_interface((1/3,2/3),self.ops['m_001'],expected_target=(.0,.0))
    def test_fractional_cartesian_distinction(self):
        # A hexagonal in-plane operation has a non-orthogonal fractional 2x2.
        found=False
        for op in self.ops.values():
            rf=np.asarray(op.matrix_fractional,float)
            rc=np.asarray(op.matrix_cartesian,float)
            if not np.allclose(rf,rc,atol=1e-8):
                found=True
                c=connect_interface((.2,.3),op)
                self.assertTrue(np.allclose(c.rotation_cartesian,rc))
                self.assertTrue(np.allclose(c.rotation_fractional,rf))
                self.assertTrue(np.allclose(rc.T@rc,np.eye(3),atol=1e-8))
                break
        self.assertTrue(found,'hex fractional and Cartesian matrices unexpectedly identical')

if __name__=='__main__':unittest.main()
