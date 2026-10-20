"""Regression gates for parity-resolved weak interlayer hybridization.

The numerical values are for synthetic 8-band model, 0K and the static
response protocol only.  Neither PRL driven DC nor physical BSF is tested.
"""
import unittest
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.channel_perturbation import evaluate_scaled_channels
from scripts.run_stage3_sign_reversal import spec_at

ROOT=Path(__file__).resolve().parents[1]

class TestStage5Physics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k,cls.w,_=valley_quadrature(SlidingBilayer().model().reciprocal,29,1.5,25,48)

    def response(self,u,t,sa=1,ca=1):
        return evaluate_scaled_channels(replace(spec_at(u),tunneling_eV=t),self.k,self.w,
                                 same_scale=sa,cross_scale=ca)

    def test_independent_new_kernel_matches_frozen_stage4(self):
        from bsf_nmee.response_mechanism import evaluate_mechanism
        spec=spec_at(.4)
        new=evaluate_scaled_channels(spec,self.k,self.w)['gamma']
        frozen=evaluate_mechanism(spec,self.k,self.w)['gamma']
        self.assertLess(abs(new-frozen),1e-12)

    def test_exact_global_layer_gauge_and_mixed_derivative(self):
        t=.005;u=.1
        pp=self.response(u,t,1,1)['gamma']
        pm=self.response(u,t,1,-1)['gamma']
        mp=self.response(u,t,-1,1)['gamma']
        mm=self.response(u,t,-1,-1)['gamma']
        self.assertLess(abs(pp-mm),1e-11)
        self.assertLess(abs(pm-mp),1e-11)
        sc4=(pp-pm-mp+mm)/(4*t*t)
        sc2=(pp-pm)/(2*t*t)
        self.assertAlmostEqual(sc4,sc2,places=7)
        self.assertLess(sc4,-.02)

    def test_weak_mixed_coefficients_negative_at_both_slidings(self):
        for u in (.1,.4):
            mixed=[]
            for t in (.0025,.005):
                a=self.response(u,t)['gamma'];b=self.response(u,t,ca=-1)['gamma']
                mixed.append((a-b)/(2*t*t))
            self.assertTrue(all(x<0 for x in mixed))
            self.assertLess(abs(mixed[0]-mixed[1]),.003)

    def test_physical_sign_reversal_is_even_odd_competition(self):
        for u,expect in ((.1,-1),(.4,1)):
            t=.020
            gp=self.response(u,t)['gamma'];gm=self.response(u,t,ca=-1)['gamma']
            even=(gp+gm)/2;odd=(gp-gm)/2
            self.assertGreater(even,0)
            self.assertLess(odd,0)
            self.assertEqual(int(np.sign(gp)),expect)
            self.assertAlmostEqual(gp,even+odd,places=14)

    def test_fixed_slide_hybridization_driven_sign_reversal_gapped(self):
        a=self.response(.2,.015);b=self.response(.2,.020)
        self.assertLess(a['gamma'],-1e-6)
        self.assertGreater(b['gamma'],1e-7)
        self.assertGreater(min(a['sampled_min_mu_distance_eV'],b['sampled_min_mu_distance_eV']),.002)

    def test_bz_normalization_and_frozen_dataset(self):
        p=ROOT/'results/sliding_stage5_weak_coupling.json'
        r=json.loads(p.read_text())
        self.assertLess(abs(r['coarse']['meta']['normalization_error'])/r['coarse']['meta']['bz_area_weight'],.002)
        self.assertLess(abs(r['fine']['meta']['normalization_error']),.002)
        self.assertTrue(r['four_quadrant']['finite_mixed_coefficient']<0)
        for u in ('0.10','0.20','0.40'):
            self.assertEqual(len(r['data'][u]['sweep']),5)
            rows=r['data'][u]['sweep']
            for row in rows:
                self.assertAlmostEqual(row['plus'],row['even_under_cross_sign']+row['odd_under_cross_sign'],places=13)

if __name__=='__main__':unittest.main()
