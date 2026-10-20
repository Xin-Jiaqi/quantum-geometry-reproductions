import unittest
import numpy as np
from igmc_lab.toolkit_bridge import (read_toolkit_snapshot,closure_order,tensor_basis,canonical_component_basis,
    classify_basis,group_invariance_max,toolkit_live_dimensions,character_dimension)

class AllMagneticPointGroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.groups=read_toolkit_snapshot()

    def test_122_mpg_exact_name_index_and_orders(self):
        self.assertEqual(len(self.groups),122)
        for group in self.groups:
            self.assertEqual(closure_order(group),group.order,group.symbol)

    def test_244_neumann_nullspaces_and_invariance(self):
        for g in self.groups:
            for kind in ('C','D'):
                b=tensor_basis(g,kind)
                self.assertLess(group_invariance_max(g,b,kind),1e-6,(g.symbol,kind))
                self.assertLessEqual(b.shape[1],9)
                self.assertIn(classify_basis(b),{'forbidden','both','longitudinal_only','transverse_only'})

    def test_244_character_average_independent_of_SVD(self):
        for g in self.groups:
            for kind in ('C','D'):
                self.assertEqual(character_dimension(g,kind),tensor_basis(g,kind).shape[1],(g.symbol,kind))
            self.assertEqual(character_dimension(g,'C',protocol='paper_printed_eq14'),
                             tensor_basis(g,'C',protocol='paper_printed_eq14').shape[1])

    def test_PT_pure_T_and_inversion_reference_cases(self):
        dims={g.symbol:tuple(tensor_basis(g,ch).shape[1] for ch in ('C','D')) for g in self.groups}
        self.assertEqual(dims['1'],(9,9))
        self.assertEqual(dims["11'"],(9,0))
        self.assertEqual(dims['-1'],(0,0))
        self.assertEqual(dims["-1'"],(0,9))
        self.assertEqual(dims["4mm1'"],(1,0))
        self.assertEqual(dims["4m'm'"],(1,2))

    def test_SMTables_I_II_selected_sector_shapes(self):
        groups={g.symbol:g for g in self.groups}
        for label in ('4mm',"4mm1'","4m'm'",'3m',"3m1'","6mm1'"):
            self.assertEqual(classify_basis(tensor_basis(groups[label],'C')),'transverse_only',label)
        self.assertEqual(classify_basis(tensor_basis(groups["4m'm'"],'D')),'longitudinal_only')
        self.assertEqual(classify_basis(tensor_basis(groups["4mm1'"],'D')),'forbidden')
        self.assertEqual(classify_basis(tensor_basis(groups["4mm"],'D')),'transverse_only')

    def test_toolkit_comparison_if_pip_installed(self):
        try:import group_theory_operations
        except ImportError:self.skipTest('external toolkit not installed; pinned generators independently checked')
        self.assertEqual(toolkit_live_dimensions(self.groups),[])

if __name__=='__main__':unittest.main()

class SourceSettingConversionTests(unittest.TestCase):
    def test_monoclinic_unique_axis_y_for_main_sm_table_I(self):
        # The toolkit's 2 operation is 2_z, while the SM Table I column uses
        # conventional monoclinic unique axis y. This is not a tensor mismatch.
        groups={g.symbol:g for g in read_toolkit_snapshot()}
        swap_yz=np.array([[1.,0,0],[0,0,1],[0,1,0]])
        b=tensor_basis(groups['2'],'C')
        masks=np.array([swap_yz@b[:,i].reshape(3,3)@swap_yz.T for i in range(b.shape[1])])
        source_allowed={(0,0),(0,2),(1,1),(2,0),(2,2)}
        for i in range(3):
            for j in range(3):
                self.assertEqual(bool(np.max(abs(masks[:,i,j]))>1e-8), (i,j) in source_allowed)

    def test_monoclinic_mirror_unique_axis_y_for_main_sm_table_I(self):
        groups={g.symbol:g for g in read_toolkit_snapshot()}
        swap_yz=np.array([[1.,0,0],[0,0,1],[0,1,0]])
        b=tensor_basis(groups['m'],'C')
        masks=np.array([swap_yz@b[:,i].reshape(3,3)@swap_yz.T for i in range(b.shape[1])])
        source_allowed={(0,1),(1,0),(1,2),(2,1)}
        for i in range(3):
            for j in range(3):
                self.assertEqual(bool(np.max(abs(masks[:,i,j]))>1e-8), (i,j) in source_allowed)

class SourceTableIIISummaryTests(unittest.TestCase):
    def test_longitudinal_hall_mixed_column_counts(self):
        from collections import Counter
        groups=read_toolkit_snapshot()
        expected={'C':{'longitudinal_only':24,'transverse_only':18,'both':16,'forbidden':64},
                  'D':{'longitudinal_only':23,'transverse_only':19,'both':16,'forbidden':64}}
        for channel in 'CD':
            actual=dict(Counter(classify_basis(tensor_basis(g,channel)) for g in groups))
            self.assertEqual(actual,expected[channel])

    def test_printed_eq14_different_C_dimension_for_58_MP_points(self):
        gs=read_toolkit_snapshot()
        inconsistent=[g for g in gs if tensor_basis(g,'C').shape[1]!=
           tensor_basis(g,'C',protocol='paper_printed_eq14').shape[1]]
        self.assertEqual(len(inconsistent),58)
        self.assertIn("4mm1'",{g.symbol for g in inconsistent})


class CanonicalTableGeneratorTests(unittest.TestCase):
    def test_all_244_pivot_basis_respects_each_group(self):
        for group in read_toolkit_snapshot():
            for kind in 'CD':
                raw=tensor_basis(group,kind)
                concrete=canonical_component_basis(raw)
                self.assertEqual(concrete.shape,raw.shape)
                self.assertLess(group_invariance_max(group,concrete,kind),1e-6)

    def test_4mprimemprime_has_expected_SM_symbolic_shape(self):
        g=next(g for g in read_toolkit_snapshot() if g.symbol=="4m'm'")
        canonical=canonical_component_basis(tensor_basis(g,'D'))
        tensors=[canonical[:,i].reshape(3,3) for i in range(2)]
        self.assertTrue(any(np.allclose(t,np.diag([1,1,0])) for t in tensors))
        self.assertTrue(any(np.allclose(t,np.diag([0,0,1])) for t in tensors))
