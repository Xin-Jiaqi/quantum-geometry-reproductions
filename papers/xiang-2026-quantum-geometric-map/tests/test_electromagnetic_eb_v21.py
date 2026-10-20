"""v21: exact finite-q mixed electric scalar / Zeeman and Peierls-gauge controls."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.electromagnetic_eb_v21 import (EBSpec,build,zeeman_operator,
 finite_zeeman_slope,liouville_slope,mixed_ward,gauge_contact,gauge_mixed,
 derivative_equilibrium_density)
from qgmap.real_frequency_continuity_v20 import ExactSystem,fermi

class TestElectromagneticEB(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.p=EBSpec(magnetic_step=2e-4)
 def test_zeeman_hermitian(self):
  z=zeeman_operator(self.p.system)
  self.assertLess(np.max(abs(z-z.conj().T)),1e-13)
 def test_charge_zeeman_commutation(self):
  _,z,r,rm,j,jm,Q,_=build(self.p)
  self.assertLess(np.max(abs(z@r-r@z)),1e-13)
 def test_current_independent_of_onsite_B(self):
  h,z,r,rm,j,jm,Q,_=build(self.p)
  self.assertLess(np.max(abs(-((h+.23*z)@r-r@(h+.23*z))/Q-j)),1e-12)
 def test_density_frechet_diagonal(self):
  h,z,*_=build(self.p);e,u=np.linalg.eigh(h)
  derivative=derivative_equilibrium_density(e,u,z,self.p.system.mu,self.p.system.temperature)
  d=np.diag(derivative)
  ze=np.diag(u.conj().T@z@u)
  f=fermi(e,self.p.system.mu,self.p.system.temperature)
  self.assertLess(np.max(abs(d+f*(1-f)*ze/self.p.system.temperature)),1e-12)
 def test_equilibrium_derivative_is_nonzero(self):
  h,z,*_=build(self.p);e,u=np.linalg.eigh(h)
  v=derivative_equilibrium_density(e,u,z,self.p.system.mu,self.p.system.temperature)
  self.assertGreater(np.linalg.norm(v),.02)
 def test_analytic_vs_fd_current_density(self):
  self.assertLess(abs(liouville_slope(self.p)-finite_zeeman_slope(self.p)),2e-7)
 def test_analytic_vs_fd_density_density(self):
  self.assertLess(abs(liouville_slope(self.p,'density','density')-finite_zeeman_slope(self.p,'density','density')),5e-7)
 def test_analytic_vs_fd_density_current(self):
  self.assertLess(abs(liouville_slope(self.p,'density','current')-finite_zeeman_slope(self.p,'density','current')),2e-7)
 def test_analytic_vs_fd_current_current(self):
  self.assertLess(abs(liouville_slope(self.p,'current','current')-finite_zeeman_slope(self.p,'current','current')),2e-7)
 def test_magnetic_fd_second_order_convergence(self):
  exact=liouville_slope(self.p)
  err_big=abs(finite_zeeman_slope(replace(self.p,magnetic_step=1e-3))-exact)
  err_small=abs(finite_zeeman_slope(replace(self.p,magnetic_step=5e-4))-exact)
  self.assertLess(err_small/err_big,.27)
 def test_density_ward(self):
  self.assertLess(mixed_ward(self.p)['density_source_ward_abs'],1e-11)
 def test_current_ward_with_contact(self):
  r=mixed_ward(self.p)
  self.assertLess(r['current_source_ward_abs'],1e-11)
  self.assertGreater(abs(r['contact_magnetic_derivative']),1e-5)
 def test_finite_field_ward(self):
  r=mixed_ward(self.p,'finite_difference')
  self.assertLess(max(r['density_source_ward_abs'],r['current_source_ward_abs']),1e-9)
 def test_equilibrium_term_is_material(self):
  full=liouville_slope(self.p)
  truncated=liouville_slope(self.p,include_equilibrium_derivative=False)
  self.assertGreater(abs(full-truncated),.1)
 def test_label_frozen_branch_is_not_equilibrium(self):
  full=finite_zeeman_slope(self.p)
  frozen=finite_zeeman_slope(self.p,occupation_protocol='frozen_label')
  self.assertGreater(abs(full-frozen),.1)
 def test_pure_gauge_static_param_contact(self):
  x=gauge_contact(self.p)
  self.assertLess(abs(x['sum']),1e-10)
  self.assertGreater(abs(x['paramagnetic']),1)
 def test_pure_gauge_finite_difference(self):
  x=gauge_contact(self.p)
  self.assertLess(abs(x['full_finite_gauge_derivative']),1e-8)
 def test_pure_gauge_mixed_order(self):
  x=gauge_mixed(self.p)
  self.assertLess(abs(x['total_mixed']),1e-9)
  self.assertLess(abs(x['finite_gauge_mixed']),2e-8)
  self.assertGreater(abs(x['paramagnetic_mixed']),1e-5)
 def test_other_momentum_direction(self):
  p=replace(self.p,direction=0)
  self.assertLess(abs(liouville_slope(p)-finite_zeeman_slope(p)),3e-7)
  self.assertLess(mixed_ward(p)['current_source_ward_abs'],1e-10)
  self.assertLess(abs(gauge_contact(p)['sum']),1e-9)
 def test_no_disorder_finite_q(self):
  p=replace(self.p,system=replace(self.p.system,disorder_std=0.))
  self.assertLess(abs(liouville_slope(p)-finite_zeeman_slope(p)),3e-7)
  self.assertLess(mixed_ward(p)['density_source_ward_abs'],1e-10)
 def test_invalid_protocol(self):
  with self.assertRaises(ValueError):finite_zeeman_slope(self.p,occupation_protocol='invalid')
 def test_invalid_frequency_regulator(self):
  with self.assertRaises(ValueError):EBSpec(delta=0.)

if __name__=='__main__':unittest.main()