"""v46 physical protocol difference, full EB two orderings, original PHE still OPEN."""
import unittest
import numpy as np
from qgmap.dirac_mixed_rta_v46 import (Spec,_blocks,quadrature,current,decompose,
                                       unprepared_two_field,audit,_b_response_unprepared)
from qgmap.dirac_phe_protocol_v45 import printed_eq11_geometry

class MixedPHEProtocolsV46(unittest.TestCase):
    def test_equilibrium_density_and_source_translate_with_B(self):
        s=Spec(nr=16,nphi=32)
        x,y,_=quadrature(s.nr,s.nphi,s.cutoff)
        for b in (.1,-.04):
            one=_blocks(s,x,y,b);two=_blocks(s,x,y+b,0)
            for i in (0,3,4):np.testing.assert_allclose(one[i],two[i],atol=4e-13)

    def test_pre_equilibrated_Hall_cross_channel_field_invariant(self):
        for eta in (.2,.065,.02):
            s=Spec(eta=eta,nr=64,nphi=128)
            for component in ('xy','yx'):
                for field in (-.05,0,.05):
                    self.assertLess(abs(current(s,field,component)),1e-5)

    def test_prepared_split_matches_independent_full_B_finite_difference(self):
        for eta in (.2,.065,.02):
            d=decompose(Spec(eta=eta,nr=72,nphi=128),h=3e-5)
            self.assertLess(abs(d['decomposition_residual']),4e-8)
            self.assertLess(abs(d['complete_split']),3e-7)
            self.assertLess(abs(d['full_finite_difference']),3e-7)

    def test_large_prepared_source_and_commutator_cancel(self):
        d=decompose(Spec(eta=.065,nr=96,nphi=144),h=3e-5)
        self.assertGreater(d['source_prepare'],.075)
        self.assertLess(d['magnetic_commutator'],-.075)
        self.assertLess(abs(d['complete_split']),1e-7)

    def test_unprepared_includes_B_then_E_and_E_then_B(self):
        r=unprepared_two_field(Spec(eta=.065,nr=96,nphi=128),step=3e-5)
        self.assertGreater(r['b_then_e'],.1)
        self.assertLess(r['e_then_b'],-.07)
        self.assertGreater(r['total'],.03)
        self.assertAlmostEqual(r['total'],r['b_then_e']+r['e_then_b'],delta=1e-12)

    def test_collision_target_changes_mixed_response(self):
        s=Spec(eta=.065,nr=96,nphi=144)
        a=decompose(s,h=3e-5)
        b=unprepared_two_field(s,step=3e-5)
        self.assertLess(abs(a['complete_split']),1e-7)
        self.assertGreater(b['total'],.03)
        self.assertAlmostEqual(b['e_then_b'],a['magnetic_commutator'],delta=3e-12)

    def test_unprepared_small_tilt_Eq8_near_not_identical(self):
        s=Spec(eta=.001,nr=180,nphi=192,temperature=.03)
        x=unprepared_two_field(s,step=2e-5)
        paper=-printed_eq11_geometry(s.beta)*np.tanh(s.mu/(2*s.temperature))
        ratio=s.eta*x['total']/paper
        self.assertLess(abs(ratio-1),.02)
        self.assertGreater(abs(ratio-1),.005)

    def test_small_eta_tau_scaling_not_equivalent_to_paper(self):
        coefficients=[]
        for eta in (.065,.02,.008):
            s=Spec(eta=eta,nr=96,nphi=128)
            coefficients.append(eta*unprepared_two_field(s,step=3e-5)['total'])
        self.assertTrue(coefficients[0]<coefficients[1]<coefficients[2])
        self.assertGreater(coefficients[-1],.0024)

    def test_full_tensor_channels_and_symmetric_dissipative_limit(self):
        a=audit(Spec(),h=3e-5)
        self.assertEqual(a['scientific_status'],'PARTIAL')
        self.assertLess(abs(a['symmetric_cross_derivative']),3e-7)
        self.assertLess(abs(a['antisymmetric_cross_derivative']),3e-7)

    def test_B_driven_density_not_equal_B_equilibrium_derivative(self):
        s=Spec(nr=16,nphi=32)
        x,y,_=quadrature(s.nr,s.nphi,s.cutoff)
        prepared=(_blocks(s,x,y,1e-5)[3]-_blocks(s,x,y,-1e-5)[3])/(2e-5)
        driven=_b_response_unprepared(s,x,y)
        self.assertGreater(np.max(abs(prepared-driven)),1e-3)

    def test_invalid_inputs_rejected(self):
        for args in ({'eta':0},{'temperature':-1},{'beta':1.5},{'nr':3}):
            with self.assertRaises(ValueError):Spec(**args)
        with self.assertRaises(ValueError):unprepared_two_field(Spec(),step=-1)
        with self.assertRaises(ValueError):decompose(Spec(),h=0)

if __name__=='__main__':unittest.main()
