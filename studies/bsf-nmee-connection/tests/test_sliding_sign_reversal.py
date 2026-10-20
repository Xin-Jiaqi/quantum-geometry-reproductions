"""Checks for stage-3 SIGN and symmetry, not PRL DC or ferroelectric stability."""
import unittest
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingBilayer,SlidingSpec
from bsf_nmee.valley_quadrature import valley_quadrature
from scripts.sliding_static_scan import response_weighted
from scripts.run_stage3_sign_reversal import spec_at,hopping_norm_ratio


class TestContinuousSliding(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.k,cls.w,cls.meta=valley_quadrature(
            SlidingBilayer().model().reciprocal,41,1.5,35,64)

    def test_sign_change_at_fixed_gap_and_high_symmetry_zero(self):
        responses=[]
        for u in (0.,.1,.2,.4):
            r=response_weighted(spec_at(u),self.k,self.w)
            self.assertGreater(r['min_mu_to_band_eV'],.001)
            responses.append(r['response']['total_x'])
        self.assertLess(abs(responses[0]),1e-11)
        self.assertLess(responses[1],-1e-6)
        self.assertGreater(responses[2],1e-7)
        self.assertGreater(responses[3],1e-5)

    def test_partner_odd_and_equal_valleys(self):
        p=response_weighted(spec_at(.4),self.k,self.w)
        n=response_weighted(spec_at(-.4),self.k,self.w)
        self.assertLess(abs(p['response']['total_x']+n['response']['total_x']),1e-10)
        self.assertLess(abs(p['response']['top_x']+n['response']['bottom_x']),1e-10)
        off=self.meta['n_exterior_kept'];nv=self.meta['n_valley_total']//2
        k1=self.k[off:off+nv];w1=self.w[off:off+nv]
        k2=self.k[off+nv:];w2=self.w[off+nv:]
        v1=response_weighted(spec_at(.4),k1,w1)['response']['total_x']
        v2=response_weighted(spec_at(.4),k2,w2)['response']['total_x']
        self.assertLess(abs(v1-v2),1e-9)
        self.assertGreater(v1,0.)

    def test_matched_average_hopping_norm_same_sublattice_is_suppressed(self):
        spec=spec_at(.4)
        n=hopping_norm_ratio(spec)
        sp_same=SlidingSpec(shift=spec.shift,tunneling_eV=spec.tunneling_eV*n['same_amplitude_rescale'],cross_sublattice=False)
        full=response_weighted(spec,self.k,self.w)
        same=response_weighted(sp_same,self.k,self.w)
        self.assertGreater(same['min_mu_to_band_eV'],.001)
        self.assertGreater(abs(full['response']['total_x']),100*abs(same['response']['total_x']))
        self.assertAlmostEqual(n['total_norm_sq_eV2'],
                               n['same_norm_sq_eV2']+n['cross_norm_sq_eV2'],places=13)

if __name__=='__main__':unittest.main()
