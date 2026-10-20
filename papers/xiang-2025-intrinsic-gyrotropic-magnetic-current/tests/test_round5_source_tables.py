"""Independent source-shape fidelity tests for visually transcribed SM Tables I/II."""
from __future__ import annotations
import unittest
import numpy as np
from igmc_lab.paper_table_source import (C_SOURCE_ROWS,D_SOURCE_ROWS,row_lookup,compare_source_table,
                            source_in_toolkit_setting,matrix_basis,Y_TO_Z)
from igmc_lab.toolkit_bridge import read_toolkit_snapshot,all_group_operations

class Round5SourceTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.groups=read_toolkit_snapshot()
        cls.entries,cls.summary=compare_source_table(cls.groups)

    def test_transcription_count_and_uniqueness(self):
        self.assertEqual(len(C_SOURCE_ROWS),11)
        self.assertEqual(len(D_SOURCE_ROWS),12)
        self.assertEqual(sum(len(r.symbols) for r in C_SOURCE_ROWS),58)
        self.assertEqual(sum(len(r.symbols) for r in D_SOURCE_ROWS),58)
        self.assertEqual(len(row_lookup()),116)

    def test_no_group_silently_skipped(self):
        self.assertEqual(len(self.entries),244)
        self.assertEqual({x['symbol'] for x in self.entries},set(g.symbol for g in self.groups))
        self.assertEqual({x['channel'] for x in self.entries},{'C','D'})

    def test_243_tensor_spaces_match_exactly(self):
        self.assertEqual(self.summary['matched_mpg_response_pairs'],243)
        bad=self.summary['mismatches']
        self.assertEqual(len(bad),1)
        self.assertEqual((bad[0]['channel'],bad[0]['symbol']),('D',"-42'm'"))
        self.assertAlmostEqual(max(x['max_projector_residual'] for x in self.entries
                                   if x['status']=='PASS'),0,places=7)

    def test_all_conduction_shapes_source_exact(self):
        self.assertTrue(all(x['status']=='PASS' for x in self.entries if x['channel']=='C'))

    def test_source_exception_is_a_missing_constraint_only(self):
        group=next(g for g in self.groups if g.symbol=="-42'm'")
        self.assertEqual(group.number,52)
        data=self.summary['mismatches'][0]
        self.assertEqual((data['sm_dimension'],data['independently_calculated_dimension']),(2,1))
        from igmc_lab.toolkit_bridge import tensor_basis
        b=tensor_basis(group,'D')
        expect=np.zeros((9,1));expect[[1,3],0]=1
        # The unit-normalized invariant form is [[0,1,0],[1,0,0],[0,0,0]].
        normalized=b[:,0]/b[1,0]
        np.testing.assert_allclose(normalized,expect[:,0],atol=1e-10)

    def test_independent_analytic_constraints(self):
        # No numerical nullspace needed: use two explicit group generators in
        # closed 2x2 tensor algebra. The S4 rotation swaps xy and yx; C2x T
        # removes xx and yy for D, leaving a symmetric single parameter.
        R_s4=np.array([[0,1,0],[-1,0,0],[0,0,-1]])
        R_c2t=np.diag([1,-1,-1])
        T=np.array([[0,1,0],[1,0,0],[0,0,0]])
        for rot,eta in ((R_s4,1),(R_c2t,-1)):
            transformed=eta*round(np.linalg.det(rot))*rot@T@rot.T
            np.testing.assert_allclose(transformed,T,atol=1e-12)
        not_allowed=np.array([[0,1,0],[0,0,0],[0,0,0]])
        self.assertFalse(np.allclose(-R_s4@not_allowed@R_s4.T,not_allowed))

    def test_monoclinic_unique_axis_transform_explicit(self):
        # Under source y -> toolkit z, xx/xz/yy/zx/zz becomes
        # xx/xy/yx/yy/zz, without affecting the dimension of the space.
        row=next(r for r in C_SOURCE_ROWS if r.id=='C02')
        b=source_in_toolkit_setting(row)
        self.assertEqual(b.shape,(9,5))
        positions={i for i in range(9) if np.linalg.norm(b[i,:])>0}
        self.assertEqual(positions,{0,1,3,4,8})
        np.testing.assert_allclose(Y_TO_Z@Y_TO_Z.T,np.eye(3))

if __name__=='__main__':unittest.main()
