"""Regression: four-band internal optical/spin-electric terms can have wrong sign."""
import unittest
from dataclasses import replace
import numpy as np
from bsf_nmee.sliding_tight_binding import SlidingSpec, SlidingBilayer
from bsf_nmee.spectral_sector import band_sector_kernel, valley_sector_integrals


class SpectralSectorAuditTests(unittest.TestCase):
    def test_bloch_velocity_derivative(self):
        model=SlidingBilayer(replace(SlidingSpec(),shift=(.1*.23,.1*.41),
                                     tunneling_eV=.005))
        k=np.array([[.13,.07]])
        delta=1.e-5
        for j in (0,1):
            basis=np.eye(2)[j]*delta
            num=(model.matrices(k+basis)-model.matrices(k-basis))/(2*delta)
            exact=model.matrices(k,axis=j)
            self.assertLess(np.max(np.abs(num-exact)),1.e-8)

    def test_four_band_only_can_have_wrong_sign(self):
        for u,t in ((.1,.02),(.4,.005)):
            spec=replace(SlidingSpec(),shift=(u*.23,u*.41),
                         tunneling_eV=t)
            row=valley_sector_integrals(spec, radius=.5,radial_n=16,angular_n=48)
            self.assertAlmostEqual(row["full"],row["middle_four_only"]+
                                   row["remote_virtual_transitions"],places=14)
            self.assertLess(row["full"]*row["middle_four_only"],0.)
            self.assertGreater(row["sampled_min_middle4_remote_gap_eV"],.002)
            self.assertGreater(row["sampled_min_mu_distance_eV"],.002)
            self.assertLess(abs(row["valley_integrals"][0][0]-
                                row["valley_integrals"][1][0]),1.e-8)

    def test_weak_t_overlaps_but_remote_is_not_negligible(self):
        u=.4
        spec=replace(SlidingSpec(),shift=(u*.23,u*.41),tunneling_eV=.005)
        row=valley_sector_integrals(spec,radius=.5,radial_n=16,angular_n=48)
        self.assertGreater(row["sampled_min_cos2_decoupled_vs_exact_middle4"],.97)
        self.assertGreater(abs(row["remote_virtual_transitions"]),abs(row["full"]))


if __name__=="__main__":
    unittest.main()
