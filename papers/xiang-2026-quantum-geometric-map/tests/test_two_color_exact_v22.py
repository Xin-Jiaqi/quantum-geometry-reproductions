"""v22 finite-q, finite-frequency two-color exact tests; not paper DC PHE."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.two_color_exact_v22 import (TwoColorSpec, frequency_channels,
 independent_sylvester_channels,static_preparation,static_eta_scan,
 time_domain_rk4,finite_amplitude_time_domain)
from qgmap.electromagnetic_eb_v21 import EBSpec,liouville_slope
from qgmap.real_frequency_continuity_v20 import ExactSystem

class TestTwoColor(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.p=TwoColorSpec()
  cls.freq=frequency_channels(cls.p)

 def test_source_order_nonzero(self):
  self.assertGreater(abs(self.freq['b_then_e']),.05)
  self.assertGreater(abs(self.freq['e_then_b']),.05)
 def test_orders_sum(self):
  self.assertLess(abs(self.freq['b_then_e']+self.freq['e_then_b']-self.freq['total']),1e-15)
 def test_independent_sylvester_b_first(self):
  self.assertLess(abs(self.freq['b_then_e']-independent_sylvester_channels(self.p)['b_then_e']),1e-11)
 def test_independent_sylvester_e_first(self):
  self.assertLess(abs(self.freq['e_then_b']-independent_sylvester_channels(self.p)['e_then_b']),1e-11)
 def test_independent_sylvester_total(self):
  self.assertLess(abs(self.freq['total']-independent_sylvester_channels(self.p)['total']),1e-11)
 def test_rho_first_order_trace(self):
  self.assertLess(abs(np.trace(self.freq['rho_e'])),1e-13)
  self.assertLess(abs(np.trace(self.freq['rho_b'])),1e-13)
 def test_second_order_trace(self):
  self.assertLess(abs(np.trace(self.freq['rho_eb'])),1e-12)
 def test_band_liouville_resolvent_residual(self):
  f=self.freq;en=f['eig']; d=en[:,None]-en[None,:]
  e=np.diag(f['occupations'].astype(complex))
  source=f['electric']@e-e@f['electric']
  x=(self.p.omega_e+1j*self.p.eta_e-d)*f['rho_e']
  self.assertLess(np.max(abs(x-source)),1e-12)
 def test_b_liouville_resolvent_residual(self):
  f=self.freq;en=f['eig']; d=en[:,None]-en[None,:]
  e=np.diag(f['occupations'].astype(complex))
  source=f['zeeman']@e-e@f['zeeman']
  x=(self.p.omega_b+1j*self.p.eta_b-d)*f['rho_b']
  self.assertLess(np.max(abs(x-source)),1e-12)
 def test_mixed_liouville_residual(self):
  f=self.freq;en=f['eig'];d=en[:,None]-en[None,:]
  x=(self.p.omega_e+self.p.omega_b+1j*(self.p.eta_e+self.p.eta_b)-d)*f['rho_eb']
  v,z=f['electric'],f['zeeman']
  y=(v@f['rho_b']-f['rho_b']@v)+(z@f['rho_e']-f['rho_e']@z)
  self.assertLess(np.max(abs(x-y)),1e-12)
 def test_time_domain_rk4(self):
  val=time_domain_rk4(self.p,t_start=-165.,dt=.04)['total']
  self.assertLess(abs(val-self.freq['total']),2e-8)
 def test_time_domain_far_past(self):
  err1=abs(time_domain_rk4(self.p,t_start=-80.,dt=.04)['total']-self.freq['total'])
  err2=abs(time_domain_rk4(self.p,t_start=-140.,dt=.04)['total']-self.freq['total'])
  self.assertLess(err2,err1/100)
 def test_exact_finite_amplitude_four_corner(self):
  ans=finite_amplitude_time_domain(self.p,.005,.0035,t_start=-155.)
  self.assertLess(abs(ans-self.freq['total']),8e-6)
 def test_finite_amplitude_quadratic_error(self):
  true=self.freq['total']
  e1=abs(finite_amplitude_time_domain(self.p,.02,.014,t_start=-155.)-true)
  e2=abs(finite_amplitude_time_domain(self.p,.01,.007,t_start=-155.)-true)
  self.assertLess(e2/e1,.27)
 def test_static_prepared_matches_v21(self):
  p=replace(self.p,omega_b=0.,omega_e=.31,eta_e=.045)
  got=static_preparation(p)['total']
  ref=liouville_slope(EBSpec(system=p.system,direction=p.direction,omega=p.omega_e,delta=p.eta_e))
  self.assertLess(abs(got-ref),1e-12)
 def test_dynamic_adiabatic_static_limit(self):
  p=replace(self.p,omega_b=0.,omega_e=.31,eta_e=.045)
  prep,rows=static_eta_scan(p,[.001,.0001,.00001,.000001])
  self.assertLess(rows[-1]['distance_to_prepared_without_population'],6e-6)
  self.assertLess(rows[-1]['distance_to_prepared_without_population'],rows[-2]['distance_to_prepared_without_population']*.12)
 def test_equilibrium_population_term_survives(self):
  p=replace(self.p,omega_b=0.,omega_e=.31,eta_e=.045)
  prep=static_preparation(p)
  self.assertGreater(abs(prep['population']),.001)
  self.assertLess(abs(prep['total']-prep['static_without_population']-prep['population']),1e-15)
 def test_omega_b_dependence(self):
  b=frequency_channels(replace(self.p,omega_b=-.17))['total']
  self.assertGreater(abs(b-self.freq['total']),.001)
 def test_eta_b_dependence(self):
  b=frequency_channels(replace(self.p,eta_b=.05))['total']
  self.assertGreater(abs(b-self.freq['total']),.001)
 def test_clean_disorder_zero(self):
  p=replace(self.p,system=replace(self.p.system,disorder_std=0.))
  a=frequency_channels(p)['total'];b=independent_sylvester_channels(p)['total']
  self.assertLess(abs(a-b),1e-10)
 def test_x_direction(self):
  p=replace(self.p,direction=0)
  a=frequency_channels(p)['total'];b=independent_sylvester_channels(p)['total']
  self.assertLess(abs(a-b),1e-10)
 def test_invalid_regulator(self):
  with self.assertRaises(ValueError):TwoColorSpec(eta_b=0)
  with self.assertRaises(ValueError):TwoColorSpec(eta_e=-.3)
 def test_invalid_direction(self):
  with self.assertRaises(ValueError):TwoColorSpec(direction=7)
 def test_static_scan_refuses_dynamic_b(self):
  with self.assertRaises(ValueError):static_eta_scan(self.p,[.1,.01])
 def test_time_domain_invalid(self):
  with self.assertRaises(ValueError):time_domain_rk4(self.p,t_start=1)
  with self.assertRaises(ValueError):finite_amplitude_time_domain(self.p,amplitude_b=0.)
 def test_static_preparation_refuses_dynamic_b(self):
  with self.assertRaises(ValueError):static_preparation(self.p)

if __name__=='__main__': unittest.main()