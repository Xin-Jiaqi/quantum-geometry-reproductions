"""Stage 6: unitary stacking symmetry forbids or allows *tensor coefficients*.

C3 and m_z are genuine operations of the synthetic microscopic model.
Tests distinguish the trace (rotational scalar) from the nematic sector.
"""
import json,unittest
from pathlib import Path
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer,lattice_geometry,M_GLOBAL
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.symmetry_scaling import static_inplane_matrix,split_parity,project_interference_tensor
from scripts.run_stage3_sign_reversal import spec_at
ROOT=Path(__file__).resolve().parents[1]

class SelectionRuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k,cls.w,_=valley_quadrature(SlidingBilayer().model().reciprocal,29,1.5,25,48)

    def test_high_symmetry_exact_mirror_channel_hamiltonians(self):
        k=np.array([[.13,.29],[8.1,-3.2],[-5.8,12.2]])
        for which in ('same','cross'):
            b=SlidingBilayer(spec_at(0))
            b.hops=[h for h in b.hops if (h[0]==h[1])==(which=='same')]
            H=b.matrices(k)
            self.assertLess(np.max(abs(M_GLOBAL@H@M_GLOBAL-H)),1e-13)

    def test_mirror_forbids_xzx_for_even_channels(self):
        mz=np.diag([1,1,-1]); ident=np.eye(3)
        T=np.zeros((3,3,3));T[0,2,0]=1.
        Q=project_interference_tensor(T,[ident,mz])
        self.assertLess(abs(Q[0,2,0]),1e-13)
        # If ONE of the two physical perturbations is odd under m_z,
        # a mixed derivative can be allowed (hypothetical irrep control).
        odd=project_interference_tensor(T,[ident,mz],characters_same=[1,-1])
        self.assertAlmostEqual(odd[0,2,0],1.)

    def test_xx_and_yy_split_mirror_and_cubic_trace(self):
        s=lambda u: replace(spec_at(u),tunneling_eV=.005)
        r01=split_parity(s(.1),self.k,self.w)
        rminus=split_parity(s(-.1),self.k,self.w)
        for key in ('full','even','mixed_odd'):
            self.assertLess(np.max(abs(r01[key]+rminus[key])),2e-10)
        # Cannot conflate the trace and the traceless sector.
        def trace(u): return np.trace(split_parity(s(u),self.k,self.w)['full'])
        def nem(u):
            M=split_parity(s(u),self.k,self.w)['mixed_odd']
            N=(M+M.T)/2-np.trace(M)*np.eye(2)/2
            return np.linalg.norm(N)
        self.assertTrue(7.3<trace(.1)/trace(.05)<8.5)
        self.assertTrue(1.8<nem(.1)/nem(.05)<2.2)

    def test_independent_C3_rotation_of_microscopic_hamiltonian(self):
        spec=replace(spec_at(.1),tunneling_eV=.005)
        a,_=lattice_geometry()
        rad=2*np.pi/3;R=np.array([[np.cos(rad),-np.sin(rad)],[np.sin(rad),np.cos(rad)]])
        shifted=(R@(np.array(spec.shift)@a))@np.linalg.inv(a)
        original=static_inplane_matrix(spec,self.k,self.w)['gamma']
        rotated=static_inplane_matrix(replace(spec,shift=tuple(shifted)),self.k,self.w)['gamma']
        self.assertLess(np.max(abs(rotated-R@original@R.T)),2e-10)

    def test_release_data_theoretical_ratios(self):
        v=json.loads((ROOT/'results/sliding_stage6_symmetry.json').read_text())
        self.assertEqual(v['tunneling_eV'],.005)
        self.assertTrue(7.3<v['empirical_scale_ratios_u01_over_u005']['trace']<8.5)
        self.assertTrue(1.8<v['empirical_scale_ratios_u01_over_u005']['mixed_nematic']<2.2)
        self.assertLess(v['rotated_shifts']['120']['matrix_residual'],2e-10)
        self.assertLess(abs(v['rotated_shifts']['60']['trace_ratio']+1),.002)
        self.assertLess(abs(v['rotated_shifts']['180']['trace_ratio']+1),1e-7)

try:
    from group_theory_operations import point_group_operations
except ImportError:
    point_group_operations=None

@unittest.skipUnless(point_group_operations is not None,'actual group-theory-operations-toolkit required')
class RealToolkitStage6Tests(unittest.TestCase):
    def test_C3_and_mirror_z_are_genuine_toolkit_operations(self):
        ops=point_group_operations('6/mmm')
        c3=[op for op in ops if np.isclose(np.asarray(op.matrix_cartesian)[2,2],1)
            and np.isclose(np.trace(np.asarray(op.matrix_cartesian)[:2,:2]),-1)]
        self.assertEqual(len(c3),2)
        mz=next(op for op in ops if op.name=='m_001')
        self.assertTrue(np.allclose(np.asarray(mz.matrix_cartesian),np.diag([1,1,-1])))
        # The TOOLKIT's fractional basis is not silently substituted for
        # this synthetic Hamiltonian's own oblique primitive basis.
        for op in c3:
            R=np.asarray(op.matrix_cartesian)
            a,_=lattice_geometry()
            induced=np.linalg.inv(a.T)@R[:2,:2]@a.T
            self.assertTrue(np.allclose(induced,np.rint(induced),atol=1e-10))
            self.assertEqual(abs(round(np.linalg.det(induced))),1)

if __name__=='__main__':unittest.main()
