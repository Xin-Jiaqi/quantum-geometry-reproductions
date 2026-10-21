"""v59: causal semi-infinite-chain lead selfenergy and independently periodic bulk Kubo checks."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.lead_selfenergy_v59 import (
 surface_g,surface_g_dx,selfenergy_and_B,ChainLead,ChainSpec,energy_transmissions,
 integrate_leads,finite_B_fd,reference_device)
from qgmap.lead_transport_v57 import LeadSpec,diamond_graph,device_response
from qgmap.embedding_magnetic_v55 import peierls_field
from qgmap.periodic_kubo_v59 import periodic_chern,qwz_h,qwz_derivatives
from qgmap.orbital_geometry_v54 import band_geometry

class SemiInfiniteLeadAndBulkKuboV59(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.spec=ChainSpec()
  cls.H,cls.pos,cls.leads=reference_device()
  cls.ref=integrate_leads(cls.H,cls.pos,cls.leads,cls.spec)

 def test_01_surface_Dyson_equation(self):
  for t in (.7,2.,7.):
   for E in (-10.,-5.,-3.,-2.,-1.,-.1,0,.1,1,2,3,5,10):
    g=surface_g(E,t,.2)
    self.assertLess(abs((E-.2-t*t*g)*g-1),1e-13)

 def test_02_retarded_causality_and_band_edges(self):
  for E in np.linspace(-10,10,121):
   g=surface_g(E,2)
   self.assertLessEqual(g.imag,1e-14)
   self.assertEqual((-2*g.imag)>1e-8,abs(E)<4)
  self.assertEqual(surface_g(4,2).imag,0)
  self.assertEqual(surface_g(-4,2).imag,0)

 def test_03_analytic_Gamma_and_real_selfenergy(self):
  for E in (-3.2,-1,0,1,3.2):
   s,ds,g,dg=selfenergy_and_B(E,0,ChainLead(0,2,.8))
   self.assertAlmostEqual(g,.8**2/4*np.sqrt(16-E*E),delta=1e-13)
   self.assertAlmostEqual(s.real,.8**2*E/8,delta=1e-13)
   self.assertEqual(ds,0j);self.assertEqual(dg,0.)

 def test_04_surface_far_energy_normalization(self):
  self.assertAlmostEqual(1e6*surface_g(1e6,2),1,delta=1e-11)
  self.assertAlmostEqual(-1e6*surface_g(-1e6,2),1,delta=1e-11)

 def test_05_selfenergy_field_derivative_finite_difference(self):
  lead=ChainLead(0,2,.7,.13,.2)
  E=.4;B=.31;h=1e-5
  s,ds,g,dg=selfenergy_and_B(E,B,lead)
  sp,_,gp,_=selfenergy_and_B(E,B+h,lead)
  sm,_,gm,_=selfenergy_and_B(E,B-h,lead)
  self.assertAlmostEqual(ds,(sp-sm)/(2*h),delta=2e-10)
  self.assertAlmostEqual(dg,(gp-gm)/(2*h),delta=2e-10)
  self.assertGreater(abs(ds),1e-3)

 def test_06_surface_energy_derivative_finite_difference(self):
  for E in (-6,-3,-1,0,.7,2.6,6):
   h=1e-5
   self.assertAlmostEqual(surface_g_dx(E,2),(surface_g(E+h,2)-surface_g(E-h,2))/(2*h),delta=5e-10)

 def test_07_full_Caroli_trace_vs_site_formula(self):
  p=peierls_field(self.H,self.pos,.31,-.87)
  E=.35;T,dT=energy_transmissions(E,p.hamiltonian,p.magnetic_vertex,self.leads,.31)
  ss=np.array([selfenergy_and_B(E,.31,l)[0] for l in self.leads])
  gg=np.array([selfenergy_and_B(E,.31,l)[2] for l in self.leads])
  G=np.linalg.inv(E*np.eye(4)-p.hamiltonian-np.diag(ss))
  for a in range(4):
   for b in range(4):
    if a==b:continue
    A=np.diag([gg[a] if i==a else 0 for i in range(4)])
    B=np.diag([gg[b] if i==b else 0 for i in range(4)])
    self.assertAlmostEqual(T[a,b],np.trace(A@G@B@G.conj().T).real,delta=1e-13)
  self.assertGreater(np.max(abs(dT)),.01)

 def test_08_optical_theorem_exact_energy_dependent_lead(self):
  E=.35;p=peierls_field(self.H,self.pos,.31,-.87)
  ss=np.array([selfenergy_and_B(E,.31,l)[0] for l in self.leads])
  gg=np.array([selfenergy_and_B(E,.31,l)[2] for l in self.leads])
  G=np.linalg.inv(E*np.eye(4)-p.hamiltonian-np.diag(ss))
  np.testing.assert_allclose(1j*(G-G.conj().T),G@np.diag(gg)@G.conj().T,atol=1e-13)

 def test_09_magnetic_gauge_every_T_and_dT(self):
  b=integrate_leads(self.H,self.pos,self.leads,self.spec,gauge='landau_y')
  np.testing.assert_allclose(self.ref['T'],b['T'],atol=2e-12)
  np.testing.assert_allclose(self.ref['dT_dB'],b['dT_dB'],atol=2e-11)
  self.assertAlmostEqual(self.ref['observables']['dHall_dB'],b['observables']['dHall_dB'],delta=1e-12)

 def test_10_full_fourprobe_five_point_field_difference(self):
  fd=finite_B_fd(self.H,self.pos,self.leads,self.spec)
  np.testing.assert_allclose(self.ref['dT_dB'],fd['dT_dB'],atol=1e-8)
  self.assertAlmostEqual(self.ref['observables']['dHall_dB'],fd['dHall_dB'],delta=1e-8)

 def test_11_B_dependent_Sigma_and_Gamma_must_be_differentiated(self):
  H,p,leads=reference_device(onsite_B_slope=.12)
  a=integrate_leads(H,p,leads)
  fd=finite_B_fd(H,p,leads)
  self.assertAlmostEqual(a['observables']['dHall_dB'],fd['dHall_dB'],delta=1e-8)
  self.assertGreater(abs(a['observables']['dHall_dB']-self.ref['observables']['dHall_dB']),.005)

 def test_12_current_conservation_and_floating_probes(self):
  self.assertLess(self.ref['current_conservation_residual'],1e-12)
  self.assertLess(self.ref['d_current_conservation_residual'],1e-12)
  self.assertLess(self.ref['observables']['max_floating_probe_current'],1e-12)

 def test_13_spinless_onsager_for_B_independent_lead_params(self):
  neg=integrate_leads(self.H,self.pos,self.leads,replace(self.spec,field=-self.spec.field))
  np.testing.assert_allclose(self.ref['T'],neg['T'].T,atol=2e-12)
  np.testing.assert_allclose(self.ref['dT_dB'],-neg['dT_dB'].T,atol=2e-11)
  self.assertAlmostEqual(self.ref['observables']['hall_ratio'],-neg['observables']['hall_ratio'],delta=1e-12)

 def test_14_finite_band_differs_but_converges_to_wideband(self):
  pos,h,sites=diamond_graph(1)
  wide=device_response(h,pos,sites,[.55]*4,LeadSpec())['observables']['hall_ratio']
  h,p,leads=reference_device(lead_hopping=2)
  low=integrate_leads(h,p,leads)['observables']['hall_ratio']
  h,p,leads=reference_device(lead_hopping=100)
  high=integrate_leads(h,p,leads)['observables']['hall_ratio']
  self.assertGreater(abs(low-wide),1e-4)
  self.assertLess(abs(high-wide),1e-5)

 def test_15_lead_outside_band_is_closed(self):
  s,ds,g,dg=selfenergy_and_B(1,0,ChainLead(0,hopping=.2,coupling=.5))
  self.assertEqual(g,0.)
  self.assertEqual(s.imag,0)
  self.assertEqual(ds,0)

 def test_16_fail_closed_band_edge_and_bad_physics(self):
  for kw in ({'hopping':0},{'coupling':0},{'site':-1},{'onsite':np.nan}):
   with self.assertRaises(ValueError):ChainLead(**{'site':0,**kw})
  with self.assertRaises(ValueError):surface_g_dx(4,2)
  with self.assertRaises(ValueError):surface_g(0,0)
  with self.assertRaises(ValueError):ChainSpec(temperature=-1)

 def test_17_independent_bulk_QWZ_charge_chern_phases(self):
  for mass,target in ((-3,0),(-1,1),(1,-1),(3,0)):
   o=periodic_chern(mass,31)
   self.assertAlmostEqual(o['chern_spectral'],target,delta=1e-7)
   self.assertAlmostEqual(o['chern_wilson'],target,delta=1e-12)

 def test_18_periodic_bulk_kmesh_convergence(self):
  coarse=periodic_chern(-1,15)
  fine=periodic_chern(-1,31)
  self.assertGreater(abs(coarse['chern_spectral']-1),abs(fine['chern_spectral']-1))
  self.assertAlmostEqual(fine['chern_wilson'],1,delta=1e-12)

 def test_19_gapless_QWZ_phase_boundary_refused(self):
  for mass in (0,-2,2):
   with self.assertRaises(ValueError):periodic_chern(mass,31)
  with self.assertRaises(ValueError):periodic_chern(-1,7)

 def test_20_independently_differentiate_bloch_QWZ(self):
  x,y=.33,.84;h=1e-5;V=qwz_derivatives(x,y)
  np.testing.assert_allclose((qwz_h(x+h,y,-1)-qwz_h(x-h,y,-1))/(2*h),V[0],atol=3e-11)
  np.testing.assert_allclose((qwz_h(x,y+h,-1)-qwz_h(x,y-h,-1))/(2*h),V[1],atol=3e-11)
  self.assertGreater(abs(band_geometry(qwz_h(x,y,-1),V).curvature[0]),.02)

 def test_21_bulk_charge_Kubo_cannot_claim_orbital_mixed_EB(self):
  self.assertIn('NOT_ORBITAL_EB',periodic_chern(-1,21)['status'])
  self.assertNotIn('B_derivative',periodic_chern(-1,21))
  self.assertIn('CONDITIONAL_DEVICE',self.ref['status'])

if __name__=='__main__':unittest.main()
