"""v20 independent exact Kubo charge conservation, thermodynamic limit distinction,
and v19 SCBA eta limit. These tests do NOT claim full original PHE acceptance."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.real_frequency_continuity_v20 import (
 ExactSystem, exact_hamiltonian, charge_operators, kinematic_ward,
 uniform_density_vs_compressibility,exact_retarded_chi,fermi)
from qgmap.limit_order_v20 import scan_limit

class TestExactDisorderRealFrequency(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.p=ExactSystem()

 def test_hermiticity(self):
  h,_,_=exact_hamiltonian(self.p)
  self.assertLess(np.max(abs(h-h.conj().T)),1e-12)

 def test_exact_lattice_commutator(self):
  h,r,_=exact_hamiltonian(self.p)
  for direction in [0,1]:
   rho,_,j,_,Q=charge_operators(h,r,direction=direction)
   self.assertLess(np.max(abs(h@rho-rho@h+Q*j)),1e-12)

 def test_current_adjoint(self):
  h,r,_=exact_hamiltonian(self.p)
  _,_,j,jbar,_=charge_operators(h,r)
  self.assertLess(np.max(abs(j.conj().T-jbar)),1e-12)

 def test_density_ward_two_directions(self):
  for d in [0,1]:
   out=kinematic_ward(self.p,direction=d)
   self.assertLess(out['ward_density_abs'],1e-11)
   self.assertLess(abs(out['contact_rho_rho']),1e-10)

 def test_current_contact_ward_two_directions(self):
  for d in [0,1]:
   out=kinematic_ward(self.p,direction=d)
   self.assertLess(out['ward_current_abs'],1e-11)
   self.assertGreater(abs(out['contact_rho_j']),5.)
   self.assertGreater(out['bare_missing_current_contact_abs'],5.)

 def test_frequency_scan(self):
  for w in [.12,.27,.50,-.23]:
   x=kinematic_ward(self.p,omega=w)
   self.assertLess(max(x['ward_density_abs'],x['ward_current_abs']),1e-10)

 def test_delta_scan(self):
  for d in [.08,.035,.01]:
   x=kinematic_ward(self.p,delta=d)
   self.assertLess(max(x['ward_density_abs'],x['ward_current_abs']),1e-10)

 def test_site_disorder_robustness(self):
  for std in [0.,.13,.35]:
   x=kinematic_ward(replace(self.p,disorder_std=std))
   self.assertLess(max(x['ward_density_abs'],x['ward_current_abs']),1e-10)

 def test_uniform_dynamic_zero(self):
  x=uniform_density_vs_compressibility(self.p)
  self.assertLess(abs(x['dynamic_uniform_density']),1e-12)

 def test_thermodynamic_compressibility_finite(self):
  x=uniform_density_vs_compressibility(self.p)
  self.assertGreater(x['equilibrium_compressibility'],0.)
  self.assertLess(abs(x['equilibrium_compressibility']-x['compressibility_finite_difference']),1e-7)

 def test_temperature_dependence(self):
  for T in [.07,.11,.18]:
   x=uniform_density_vs_compressibility(replace(self.p,temperature=T))
   self.assertGreater(x['equilibrium_compressibility'],0.)
   self.assertLess(abs(x['equilibrium_compressibility']-x['compressibility_finite_difference']),1e-6)

 def test_invalid_direction(self):
  with self.assertRaises(ValueError):kinematic_ward(self.p,direction=3)

 def test_invalid_delta(self):
  with self.assertRaises(ValueError):kinematic_ward(self.p,delta=0.)

 def test_invalid_zero_momentum_current(self):
  h,r,_=exact_hamiltonian(self.p)
  with self.assertRaises(ValueError):charge_operators(h,r,q=0.)

class TestBornLimit(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.r=scan_limit(n=41,eta_values=(.02,.005,.001))

 def test_finite_probe_regulator_goes_to_zero(self):
  x=self.r['values']
  a=[v['fixed_probe_regulator_leak_abs'] for v in x]
  self.assertGreater(a[0],a[1]);self.assertGreater(a[1],a[2])

 def test_uniform_static_nonzero_leak(self):
  x=self.r['values']
  self.assertGreater(x[-1]['uniform_static_regulator_leak_abs'],.01)
  self.assertGreater(x[-1]['uniform_static_gamma_norm'],10.)

 def test_conserving_matrix_identity_and_finite_scattering(self):
  for v in self.r['values']:
   self.assertLess(v['fixed_probe_ward_abs'],2e-9)
   self.assertLess(v['uniform_static_ward_abs'],2e-9)
   self.assertGreater(v['disorder_scba_broadening_min'],.001)

 def test_prohibits_zero_probe_as_finite_scan(self):
  with self.assertRaises(ValueError):scan_limit(q=0.,omega=0.)

if __name__=='__main__': unittest.main()
