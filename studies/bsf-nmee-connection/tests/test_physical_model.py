"""End-to-end test: actual toolkit connection vs independent 8-band response."""
import sys
import unittest
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(BASE),str(BASE/'scripts')]
from bsf_nmee import connect_interface, transform_layer_responses
from mirror_covariance_nmee import calculate

try:
    from group_theory_operations import point_group_operations
except ImportError:
    point_group_operations=None

@unittest.skipUnless(point_group_operations is not None,'Requires actual upstream group-theory toolkit')
class PhysicalHamiltonianCovarianceTests(unittest.TestCase):
    def test_connection_predicts_independently_calculated_8band_response(self):
        mirror=next(op for op in point_group_operations('6/mmm') if op.name=='m_001')
        bridge=connect_interface((1/3,2/3),mirror,expected_target=(2/3,1/3))
        dataset=calculate(mesh_n=11)
        for branch in ('printed','full_static'):
            src={}
            for layer in ('top','bottom'):
                matrix=np.zeros((3,3,3))
                matrix[0,2,0]=dataset['results']['tau1'][layer+'_x'][branch]
                src[layer]=matrix
            prediction=transform_layer_responses(src,bridge)
            for layer in ('top','bottom'):
                obtained=dataset['results']['tau2'][layer+'_x'][branch]
                self.assertLess(abs(obtained-prediction[layer][0,2,0]),2e-11)

if __name__=='__main__':unittest.main()
