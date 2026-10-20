"""Scientific regression tests for the independently generated 8-band study."""
import unittest
from dataclasses import replace
from math import pi
from bsf_igmc_micro.model import (Parameters, hamiltonian, connection_error,
    geometric_response, layer_polarization, wilson_zeeman_slope, overlaps)
from bsf_igmc_micro.symmetry import screen
import numpy as np

P=Parameters(mass=.85)

class TestPhysicalModel(unittest.TestCase):
  def test_hilbert_space_and_hermitian(self):
    h,v=hamiltonian(.34,.51,.19,p=P)
    self.assertEqual(h.shape,(8,8))
    self.assertTrue(np.allclose(h,h.conj().T,atol=1e-13))
    self.assertEqual(set(v),{'x','y'})

  def test_geometric_hopping_is_periodic(self):
    T,_,_=overlaps(.37,-.9,.19,P)
    Q,_,_=overlaps(.37+2*pi,-.9-2*pi,.19,P)
    self.assertLess(np.max(np.abs(T-Q)),1e-13)

  def test_no_phenomenological_delta_allowed(self):
    with self.assertRaisesRegex(ValueError,'layer-odd onsite'):
      hamiltonian(.2,.4,.19,p=P,onsite_layer=.1)

  def test_both_partner_connections_and_mxprime(self):
    for b in ('AF','FM'):
      errs=connection_error(.19,b,P,N=3)
      self.assertLess(errs['partner'],1e-12)
      self.assertLess(errs['internal_mxT'],1e-12)
      self.assertLess(errs['hermiticity'],1e-12)
      self.assertLess(errs['derivative_x'],1e-8)
      self.assertLess(errs['derivative_y'],1e-8)

  def test_toolkit_pinned_character(self):
    x=screen()
    self.assertTrue(x['internal_state_groups']['mxprime']['D']['xx_allowed'])
    self.assertEqual(x['connection_pairs']['mz']['Dxx_switch_sign'],-1)
    self.assertEqual(x['connection_pairs']['mzprime']['Dxx_switch_sign'],+1)

  def test_half_filled_insulator_and_partner(self):
    for b,mu,sign in [('AF',0.,-1),('FM',-.03,+1)]:
      a=geometric_response(+.19,b,P,N=15,mu=mu)
      c=geometric_response(-.19,b,P,N=15,mu=mu)
      self.assertGreater(a['indirect_gap'],.06)
      self.assertEqual(a['occ_range'],[4,4])
      self.assertGreater(abs(a['Dxx_over_omega']),.001)
      self.assertAlmostEqual(a['Dxx_over_omega'],sign*c['Dxx_over_omega'],delta=3e-12)
      self.assertLess(abs(a['Cxx']),1e-10)

  def test_emergent_layer_dipole_and_switching(self):
    for b,mu in [('AF',0.),('FM',-.03)]:
      a=layer_polarization(.19,b,P,N=11,mu=mu)
      c=layer_polarization(-.19,b,P,N=11,mu=mu)
      self.assertGreater(abs(a),.05)
      self.assertAlmostEqual(a,-c,delta=1e-11)

  def test_nonsliding_state_zero_af(self):
    r=geometric_response(0.,'AF',P,N=13,mu=0.)
    self.assertLess(abs(r['Dxx_over_omega']),1e-11)
    self.assertLess(abs(layer_polarization(0.,'AF',P,N=13)),1e-11)

  def test_spin_conserving_interlayer_ablation(self):
    r0=geometric_response(.19,'AF',replace(P,inter_spin=0.),N=17,mu=0.)
    r1=geometric_response(.19,'AF',P,N=17,mu=0.)
    self.assertGreater(abs(r1['Dxx_over_omega']),20*abs(r0['Dxx_over_omega']))
    self.assertGreater(r0['indirect_gap'],.1)

  def test_wilson_independence(self):
    for b,mu in [('AF',0.),('FM',-.03)]:
      # Wilson uses an independently diagonalized finite-B polarization,
      # not the IGMC band's interband pair-kernel implementation.
      a=geometric_response(.19,b,P,N=17,mu=mu)['Dxx_over_omega']
      w=wilson_zeeman_slope(.19,b,P,N=17,field=3e-4)
      self.assertLess(abs(a-w),.0002)

  def test_mesh_convergence(self):
    for b,mu in [('AF',0.),('FM',-.03)]:
      a=geometric_response(.19,b,P,N=17,mu=mu)['Dxx_over_omega']
      d=geometric_response(.19,b,P,N=25,mu=mu)['Dxx_over_omega']
      self.assertLess(abs(a-d),1e-6)

  def test_antiunitary_mirror_defines_only_equivalence(self):
    a=geometric_response(0.,'FM',P,N=11,mu=-.03)
    self.assertGreater(abs(a['Dxx_over_omega']),.001)
    # Being mirror-connected by mz*T does not require zero D at u=0.

if __name__=='__main__':unittest.main()
