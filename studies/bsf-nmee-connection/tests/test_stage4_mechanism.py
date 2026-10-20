"""Model-level physics gates, not an assertion of universal NMEE physics."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer
from bsf_nmee.valley_quadrature import valley_quadrature
from bsf_nmee.response_mechanism import evaluate_mechanism
from scripts.run_stage3_sign_reversal import spec_at
from scripts.sliding_static_scan import response_weighted

class TestOrbitalChannelPhysics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k,cls.w,_=valley_quadrature(SlidingBilayer().model().reciprocal,29,1.5,25,48)

    def test_spectral_identity_vs_independent_prior_kernel(self):
        s=spec_at(.4)
        actual=evaluate_mechanism(s,self.k,self.w)
        old=response_weighted(s,self.k,self.w)
        self.assertLess(abs(actual['gamma']-old['response']['total_x']),1e-12)
        self.assertLess(abs(actual['decomposition_closure_error']),1e-13)

    def test_channel_interference_and_relative_sign_are_physical_controls(self):
        for u in (.1,.4):
            s=spec_at(u)
            full=evaluate_mechanism(s,self.k,self.w)
            cross=evaluate_mechanism(s,self.k,self.w,hopping_channel='cross_only')
            same=evaluate_mechanism(s,self.k,self.w,hopping_channel='same_only')
            reversed_cross=evaluate_mechanism(s,self.k,self.w,cross_amplitude_sign=-1)
            self.assertGreater(min(x['sampled_min_mu_distance_eV'] for x in (full,cross,same,reversed_cross)),.001)
            self.assertGreater(cross['gamma'],0)
            self.assertLess(abs(same['gamma']),1e-7)
            if u==.1:
                self.assertLess(full['gamma'],-1e-6)
                self.assertGreater(reversed_cross['gamma'],1e-6)
                self.assertLess(full['gamma']-cross['gamma']-same['gamma'],-1e-6)
            else:
                self.assertGreater(full['gamma'],2e-5)
                self.assertGreater(full['gamma']-cross['gamma']-same['gamma'],1e-6)

    def test_layer_gauge_even_t_and_weak_t_squared(self):
        s=spec_at(.1)
        t1=.0025;t2=.005
        p1=evaluate_mechanism(replace(s,tunneling_eV=t1),self.k,self.w)['gamma']
        n1=evaluate_mechanism(replace(s,tunneling_eV=-t1),self.k,self.w)['gamma']
        p2=evaluate_mechanism(replace(s,tunneling_eV=t2),self.k,self.w)['gamma']
        self.assertLess(abs(p1-n1),1e-11)
        self.assertLess(abs(p2/p1-4.),.3)

if __name__=='__main__':unittest.main()
