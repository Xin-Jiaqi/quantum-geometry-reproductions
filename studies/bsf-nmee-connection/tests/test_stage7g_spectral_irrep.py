"""Stage7g physical regression: radial cancellation hidden by total valley response."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec
from bsf_nmee.spectral_irrep_radius import (
    spectral_tensor, channel_irrep_sectors)
from bsf_nmee.spectral_sector import valley_sector_integrals


def spec_for(u,t=.0003125):
    return replace(SlidingSpec(),shift=(u*.23,u*.41),tunneling_eV=t)


class SpectralIrrepRadiusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec=spec_for(.4)
        cls.inner=channel_irrep_sectors(spec,radius=.1,nr=18,nt=48)
        cls.outer=channel_irrep_sectors(spec,radius=.75,nr=18,nt=48)

    def test_scalar_b_sign_is_remote_mediated_and_radius_dependent(self):
        a=self.inner['channels']['B']; b=self.outer['channels']['B']
        innerLL=a['per_sector']['low_to_low']['scalar']
        outerLL=b['per_sector']['low_to_low']['scalar']
        self.assertGreater(a['total']['scalar'],0)
        self.assertGreater(b['total']['scalar'],0)
        self.assertGreater(innerLL,0)
        self.assertLess(outerLL,0)
        self.assertLess(a['total']['scalar']-innerLL,0)
        self.assertGreater(b['total']['scalar']-outerLL,0)
        self.assertLess(abs(a['total']['scalar']-b['total']['scalar']) /
                        b['total']['scalar'],.02)

    def test_nematic_c_not_equivalent_to_scalar_b(self):
        b=self.outer['channels']['B']
        c=self.outer['channels']['C']
        self.assertGreater(b['total']['scalar'],0)
        self.assertLess(b['per_sector']['low_to_low']['scalar'],0)
        self.assertLess(c['total']['nematic_xx'],0)
        self.assertLess(c['per_sector']['low_to_low']['nematic_xx'],0)
        self.assertGreater(self.outer['min_mu_gap_eV'],.002)
        self.assertGreater(self.outer['min_remote_group_gap_eV'],.008)

    def test_tensor_xx_sectors_match_frozen_stage7b_engine(self):
        s=spec_for(.1,.005)
        new=spectral_tensor(s,radius=.5,nr=12,nt=36)
        old=valley_sector_integrals(s,radius=.5,radial_n=12,angular_n=36)
        np.testing.assert_allclose(new['sectors'][:,0,0],
                                   old['sum_integrals'],rtol=1e-9,atol=1e-12)

    def test_each_group_is_mirror_odd(self):
        plus=channel_irrep_sectors(spec_for(.1),radius=.75,nr=10,nt=28)
        minus=channel_irrep_sectors(spec_for(-.1),radius=.75,nr=10,nt=28)
        for ch,sector in [('B','scalar'),('C','nematic_xx')]:
            for category in plus['categories']:
                value=plus['channels'][ch]['per_sector'][category][sector]
                mirrored=minus['channels'][ch]['per_sector'][category][sector]
                self.assertLess(abs(value+mirrored),max(2e-8,1e-3*abs(value)))


if __name__=='__main__':
    unittest.main()
