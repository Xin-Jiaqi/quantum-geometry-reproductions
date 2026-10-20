"""Stage-3 false-ferroelectric-claim rejection checks."""
import unittest
from dataclasses import replace
from bsf_igmc_micro.model import Parameters, band_energy, geometric_response, layer_polarization

P=replace(Parameters(mass=.85),shell=3)

class TestElectronicSlidingGate(unittest.TestCase):
  def test_band_sum_selects_translation_identified_half_shift(self):
    # This is the frozen specific model; these inequalities are NOT a DFT claim.
    for b in ('AF','FM'):
      es=[band_energy(u,b,P,N=11) for u in (0.,.19,.5,.81,1.)]
      self.assertGreater(es[0],es[1])
      self.assertGreater(es[1],es[2])
      self.assertGreater(es[3],es[2])
      self.assertLess(abs(es[0]-es[4]),1e-8)

  def test_high_symmetry_configuration_is_nonpolar_but_not_zero_igmc_FM(self):
    a=layer_polarization(.5,'FM',P,N=11,mu=-.03)
    r=geometric_response(.5,'FM',P,N=11,mu=-.03)
    self.assertLess(abs(a),1e-7)
    self.assertGreater(r['indirect_gap'],0.15)
    self.assertGreater(abs(r['Dxx_over_omega']),.004)

  def test_high_symmetry_configuration_is_nonpolar_and_zero_igmc_AF(self):
    a=layer_polarization(.5,'AF',P,N=11,mu=0.)
    r=geometric_response(.5,'AF',P,N=11,mu=0.)
    self.assertLess(abs(a),1e-7)
    self.assertGreater(r['indirect_gap'],0.2)
    self.assertLess(abs(r['Dxx_over_omega']),1e-7)

if __name__=='__main__':unittest.main()
