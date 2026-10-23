"""v75 physically honest theory gates: no claim of entire Xiang source closure."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.elastic_spin_collision_v75 import (
 Spec,collision_data,steady,eig_physics,sign_invariance,
 independent_weighted_spectral,time_domain_independent)
from qgmap.selected_electric_pole_v74 import Spec as PS, literal_eq20,printed_eq31

class ElasticSpinCollisionV75(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.scalar=Spec(beta=0,spin_overlap=0)
  cls.chiral=Spec(beta=0,spin_overlap=1)
  cls.tilt=Spec(beta=.125,spin_overlap=1)
  cls.strong=Spec(beta=.5,spin_overlap=1)
  cls.ds={k:collision_data(s) for k,s in [('scalar',cls.scalar),('chiral',cls.chiral),('tilt',cls.tilt),('strong',cls.strong)]}
  cls.ss={k:steady(s) for k,s in [('scalar',cls.scalar),('chiral',cls.chiral),('tilt',cls.tilt),('strong',cls.strong)]}
  cls.ep={k:eig_physics(s) for k,s in [('scalar',cls.scalar),('chiral',cls.chiral),('tilt',cls.tilt),('strong',cls.strong)]}
  rng=np.random.default_rng(734)
  a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
  b=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
  cls.r=(a+a.conj().T)/2;cls.z=(b+b.conj().T)/2
  cls.e=np.array([-.8,-.21,.51,1.2]);cls.df=np.array([.3,-.11,.19,-.25])
  cls.sign=sign_invariance(cls.e,cls.r,cls.z,cls.df)

 def test_01_valid_normalized_dos(self):
  for d in self.ds.values():self.assertAlmostEqual(sum(d['weights']),1.,places=14)
 def test_02_contour_positive(self):
  for d in self.ds.values():self.assertTrue(np.all(d['radius']>0))
 def test_03_isoenergy_exact(self):
  for e in self.ep.values():self.assertLess(e['energy_shell_error'],1e-13)
 def test_04_collision_rates_positive(self):
  for d in self.ds.values():self.assertGreaterEqual(np.min(d['collision']),0.)
 def test_05_collision_detailed_balance(self):
  for e in self.ep.values():self.assertLess(e['detailed_balance_error'],1e-16)
 def test_06_number_conservation(self):
  for e in self.ep.values():self.assertLess(e['number_conservation_error'],1e-15)
 def test_07_uniform_mode_zero(self):
  for e in self.ep.values():self.assertLess(e['constant_mode_error'],1e-14)
 def test_08_unique_conserved_mode(self):
  for e in self.ep.values():
   self.assertLess(abs(e['zero_mode']),1e-12)
   self.assertGreater(e['spectral_gap'],.01)
 def test_09_positive_weighted_spectrum(self):
  for e in self.ep.values():self.assertGreaterEqual(e['spectral_gap'],0.)
 def test_10_electric_source_has_zero_density_component(self):
  for z in self.ss.values():self.assertLess(abs(z['source_weighted_mean']),1e-14)
 def test_11_response_charge_conservation(self):
  for z in self.ss.values():self.assertLess(abs(z['weighted_number']),1e-12)
 def test_12_inversion_linear_residual(self):
  for z in self.ss.values():self.assertLess(z['linear_solver_residual'],1e-12)
 def test_13_scalar_unoverlapped_tau(self):
  self.assertAlmostEqual(self.scalar.gamma*self.ss['scalar']['tau_current'].real,1,places=12)
 def test_14_helicity_doubles_transport_lifetime(self):
  self.assertAlmostEqual(self.chiral.gamma*self.ss['chiral']['tau_current'].real,2,places=12)
 def test_15_transport_rate_half_quantum_rate_in_untilted_dirac(self):
  self.assertAlmostEqual(self.ep['chiral']['spectral_gap'],self.chiral.gamma/2,places=12)
 def test_16_scalar_transport_gap(self):
  self.assertAlmostEqual(self.ep['scalar']['spectral_gap'],self.scalar.gamma,places=12)
 def test_17_tilt_breaks_universal_scalar_tau(self):
  self.assertGreater(abs(self.tilt.gamma*self.ss['tilt']['tau_current'].real-2),.02)
 def test_18_strong_tilt_reshapes_tau(self):
  self.assertLess(self.strong.gamma*self.ss['strong']['tau_current'].real,1.6)
 def test_19_nonzero_tilt_quantum_rate_angle_dependent(self):
  self.assertGreater(np.ptp(self.ds['tilt']['row_rates']),.001)
 def test_20_untilted_quantum_rate_constant(self):
  self.assertLess(np.ptp(self.ds['chiral']['row_rates']),1e-13)
 def test_21_equal_radius_is_not_elastic_when_tilted(self):
  self.assertGreater(self.ep['tilt']['equal_radius_energy_spread'],.07)
 def test_22_equal_radius_is_elastic_when_untilted(self):
  self.assertEqual(self.ep['chiral']['equal_radius_energy_spread'],0.)
 def test_23_underlying_spin_overlap_symmetric(self):
  d=self.ds['tilt'];W=d['collision'];w=d['weights'];
  self.assertLess(np.max(abs(w[:,None]*W-w[None,:]*W.T)),1e-15)
 def test_24_spectral_independent_untilted(self):
  z=independent_weighted_spectral(self.chiral)
  self.assertLess(abs(z['tau_current']-self.ss['chiral']['tau_current']),1e-11)
 def test_25_spectral_independent_tilt(self):
  z=independent_weighted_spectral(self.tilt)
  self.assertLess(abs(z['tau_current']-self.ss['tilt']['tau_current']),1e-10)
 def test_26_spectral_independent_strong(self):
  z=independent_weighted_spectral(self.strong)
  self.assertLess(abs(z['tau_current']-self.ss['strong']['tau_current']),1e-10)
 def test_27_one_zero_mode_spectral(self):
  self.assertEqual(independent_weighted_spectral(self.tilt)['zero_mode_count'],1)
 def test_28_AC_matches_independent_spectrum(self):
  for omega in (.01,.07,.27,-.12):
   a=steady(self.tilt,omega)['tau_current'];b=independent_weighted_spectral(self.tilt,omega)['tau_current']
   self.assertLess(abs(a-b),1e-10)
 def test_29_AC_reality(self):
  a=steady(self.tilt,.1)['tau_current'];b=steady(self.tilt,-.1)['tau_current']
  self.assertLess(abs(a-b.conjugate()),1e-12)
 def test_30_frequency_decreases_longitudinal_response(self):
  self.assertLess(abs(steady(self.chiral,.12)['tau_current']),abs(self.ss['chiral']['tau_current']))
 def test_31_time_domain_converges(self):
  s=Spec(beta=.125,n=24);a=time_domain_independent(s,duration=500)
  b=steady(s)
  self.assertLess(np.max(abs(a['distribution']-b['distribution'])),2e-5)
 def test_32_time_domain_conserves(self):
  z=time_domain_independent(Spec(n=24),duration=50)
  self.assertLess(abs(z['weighted_number']),1e-12)
 def test_33_source_sign_direct_vs_printed(self):
  a=self.sign['literal'];self.assertLess(abs(a['direct_eq20']+a['printed_eq31']),1e-13)
 def test_34_source_sign_is_nontrivial(self):
  self.assertGreater(abs(self.sign['literal']['direct_eq20']),.15)
 def test_35_spin_sign_flip_changes_both(self):
  self.assertLess(max(abs(x) for x in self.sign['spin_flip_both']),1e-13)
 def test_36_position_sign_flip_changes_both(self):
  self.assertLess(max(abs(x) for x in self.sign['position_flip_both']),1e-13)
 def test_37_relative_sign_survives_common_flip(self):
  self.assertLess(abs(self.sign['relative_sum']),1e-13)
 def test_38_literal_eq20_dirac_still_opposite_eq31(self):
  s=PS(nr=40,nphi=64)
  self.assertLess(abs(literal_eq20(s)+printed_eq31(s)),1e-14)
  self.assertGreater(abs(literal_eq20(s)),.002)
 def test_39_cannot_relabel_tau_as_quantum_lifetime(self):
  self.assertGreater(self.ss['chiral']['tau_current'].real-1/self.chiral.gamma,5.)
 def test_40_positivity_and_strict_input_validation(self):
  for k in ({'beta':1.},{'gamma':0.},{'energy':-1.},{'n':13},{'spin_overlap':1.2},{'beta':np.inf}):
   with self.assertRaises(ValueError):Spec(**k)
 def test_41_time_domain_bad_duration(self):
  with self.assertRaises(ValueError):time_domain_independent(Spec(n=24),0.)
 def test_42_frequency_nonfinite(self):
  with self.assertRaises(ValueError):steady(self.tilt,np.inf)
  with self.assertRaises(ValueError):independent_weighted_spectral(self.tilt,np.nan)
 def test_43_anisotropy_is_not_single_angle_independent_rate(self):
  d=self.ds['strong'];self.assertGreater(np.std(d['row_rates'])/np.mean(d['row_rates']),.1)
 def test_44_physics_scope_is_explicit(self):
  self.assertIn('NOT_FULL_PHE',self.ss['chiral']['status'])
  self.assertIn('NOT_DISORDER_FIT',self.ep['tilt']['status'])
 def test_45_tau_gamma_scaling(self):
  s1=Spec(beta=.125,gamma=.07);s2=replace(s1,gamma=.14)
  self.assertLess(abs(steady(s1)['tau_current']/2-steady(s2)['tau_current']),1e-11)
 def test_46_overlap_strength_changes_transport(self):
  s0=replace(self.tilt,spin_overlap=0.)
  self.assertGreater(abs(steady(s0)['tau_current']-self.ss['tilt']['tau_current']),8.)

if __name__=='__main__':unittest.main()
