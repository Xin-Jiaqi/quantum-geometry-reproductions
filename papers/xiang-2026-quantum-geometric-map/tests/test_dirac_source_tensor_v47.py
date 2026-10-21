"""v47: fail closed on the v46 wrong y-electric source, test real PHE symmetry."""
import unittest
import numpy as np
from qgmap.dirac_mixed_rta_v46 import Spec,quadrature,_blocks,decompose,unprepared_two_field
from qgmap.dirac_source_tensor_v47 import (
    derivative_source_check,unprepared_yx,tensor,
    independent_mixed_B_derivative,independent_finite_B_current)
from qgmap.dirac_phe_protocol_v45 import angular_eq8_geometry


class V47OriginalPaperTensor(unittest.TestCase):
    def test_source_finite_difference_catches_v46_sign_error(self):
        for beta in (.025,.125,.6):
            d=derivative_source_check(Spec(beta=beta,nr=18,nphi=32),step=2e-5)
            self.assertLess(d['max_source_y_fd_error'],2e-7)
            self.assertGreater(d['max_source_y_fd_value'],.01)

    def test_source_H_B_shift_comes_from_same_full_matrix(self):
        s=Spec(nr=16,nphi=32)
        x,y,_=quadrature(s.nr,s.nphi,s.cutoff)
        for b in (-.1,.08):
            q=_blocks(s,x,y,b);r=_blocks(s,x,y+b,0)
            np.testing.assert_allclose(q[4],r[4],atol=1e-11)

    def test_corrected_v46_prepared_terms_double_and_cancel(self):
        q=decompose(Spec(eta=.065,nr=96,nphi=144),h=3e-5)
        self.assertGreater(q['source_prepare'],.15)
        self.assertLess(q['magnetic_commutator'],-.15)
        self.assertLess(abs(q['complete_split']),2e-8)
        self.assertLess(abs(q['decomposition_residual']),2e-8)

    def test_stale_frozen_xy_total_must_be_retracted(self):
        s=Spec(eta=.065,nr=96,nphi=144)
        xy=unprepared_two_field(s,step=3e-5)
        self.assertLess(xy['total'],-.04)
        self.assertAlmostEqual(xy['total'],-.0421997868191,delta=3e-7)
        self.assertGreater(abs(xy['total']-.03377500634446),.07)

    def test_yx_channel_and_symmetry_are_not_xy_only(self):
        t=tensor(Spec(eta=.065,nr=96,nphi=144),step=3e-5)
        self.assertAlmostEqual(t['xy']['total'],-.0421997868191,delta=2e-7)
        self.assertAlmostEqual(t['yx']['total'],-.0389764539560,delta=2e-7)
        self.assertAlmostEqual(t['symmetric'],-.0405881203875,delta=3e-7)
        self.assertLess(abs(t['antisymmetric']),.002)
        self.assertEqual(t['scientific_status'],'PARTIAL')

    def test_low_broadening_symmetric_tensor_matches_Eq11_across_tilt(self):
        for beta in (.025,.08,.125,.3,.6):
            s=Spec(beta=beta,eta=.002,temperature=.03,nr=160,nphi=192)
            t=tensor(s,step=2e-5)
            self.assertLess(abs(t['ratio_to_main_eq11']-1.),2e-4)
            self.assertLess(abs(t['eta_times_antisymmetric']),4e-7)
            self.assertLess(t['eta_times_symmetric'],0)

    def test_Fermi_source_vs_paper_printed_Eq8_sign_stays_OPEN(self):
        s=Spec(eta=.002,temperature=.03,nr=160,nphi=192)
        t=tensor(s,step=2e-5)
        original_eq8,_=angular_eq8_geometry(s.beta,1)
        self.assertGreater(original_eq8,0)
        self.assertLess(t['printed_eq11_geometry'],0)
        self.assertEqual(t['scientific_status'],'PARTIAL')

    def test_independent_full_field_superoperator_matches_both_orderings(self):
        for eta in (.065,.02):
            s=Spec(eta=eta,nr=20,nphi=32)
            t=tensor(s,step=1e-4)
            for axis,key in (('y','xy'),('x','yx')):
                independent=independent_mixed_B_derivative(s,axis=axis,db=1e-4,dk=1e-4)
                self.assertAlmostEqual(independent,t[key]['total'],delta=2e-6)

    def test_direct_finite_B_current_is_not_individually_zero_in_frozen_bath(self):
        s=Spec(eta=.065,nr=20,nphi=32)
        for axis in ('x','y'):
            a=independent_finite_B_current(s,b=1e-4,axis=axis)
            b=independent_finite_B_current(s,b=-1e-4,axis=axis)
            self.assertGreater(abs(a-b),5e-6)

    def test_invalid_model_and_steps_fail_closed(self):
        with self.assertRaises(ValueError):derivative_source_check(step=-1)
        with self.assertRaises(ValueError):unprepared_yx(step=0)
        with self.assertRaises(ValueError):tensor(Spec(beta=-.2))
        with self.assertRaises(ValueError):independent_finite_B_current(Spec(),axis='z')
        with self.assertRaises(ValueError):independent_mixed_B_derivative(db=0)

if __name__=='__main__':unittest.main()
