"""Physics gates for Stage 7c: channel parity, virtual remote states, gap."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.channel_spectral_mechanism import (
    channel_mechanism, disk_channel_sectors)


def spec_at(u, t=.005):
    return replace(SlidingSpec(), shift=(.23*u,.41*u), tunneling_eV=t)


class Stage7CChannelMechanismTests(unittest.TestCase):
    def test_negative_interference_positive_cross_channel(self):
        for u in (.1,.4):
            res=channel_mechanism(spec_at(u),radius=.75,radial_n=14,angular_n=48)
            self.assertLess(res["mixed_total"],0)
            self.assertGreater(res["even_total"],0)
            self.assertGreater(res["cross_only_total"],0)
            self.assertLess(abs(res["same_only_total"]),
                            .02*abs(res["cross_only_total"]))
            self.assertGreater(res["min_mu_distance_eV"],.002)
            self.assertAlmostEqual(
                res["full"]/res["t_eV"]**2,
                res["mixed_total"]+res["even_total"],places=9)

    def test_remote_sector_decides_high_sliding_sign(self):
        res=channel_mechanism(spec_at(.4),radius=.75,radial_n=14,angular_n=48)
        self.assertLess(res["full_low_low_over_t2"],0)
        self.assertGreater(res["full_remote_over_t2"],0)
        self.assertGreater(res["full_remote_over_t2"],
                           abs(res["full_low_low_over_t2"]))
        self.assertGreater(res["full"],0)

    def test_sliding_reverse_is_odd_per_sector(self):
        kw=dict(radius=.75,radial_n=10,angular_n=32,cross_scale=-1.)
        plus=disk_channel_sectors(spec_at(.1),**kw)
        minus=disk_channel_sectors(spec_at(-.1),**kw)
        self.assertTrue(np.allclose(
            plus["sectors"],-np.asarray(minus["sectors"]),
            atol=1e-10,rtol=1e-5))

    def test_sign_reversal_bracket(self):
        kw=dict(radius=.75,radial_n=14,angular_n=48)
        below=disk_channel_sectors(spec_at(.28),**kw)
        above=disk_channel_sectors(spec_at(.31),**kw)
        self.assertLess(below["total"],0)
        self.assertGreater(above["total"],0)


if __name__=="__main__":
    unittest.main()
