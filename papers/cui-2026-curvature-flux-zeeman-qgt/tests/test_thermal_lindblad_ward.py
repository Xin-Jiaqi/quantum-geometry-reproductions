"""Round29 independent canonical/KMS, local Ward and bath gauge-vertex checks."""
import unittest
import numpy as np
from flux_lab.thermal_lindblad_ward import (
 two_site_soc,canonical_gibbs,local_lindblad,thermal_transition_jumps,
 dissipator,liouvillian,master_rhs,population_rhs,bond_current,
 continuity_terms,gauge_rotate,bath_gauge_vertex,evolved_rk4,evolved_dense)

H,P,A=two_site_soc()
BETA=.8
GAMMA=.16
G=np.diag([.47,.47,-.33,-.33]).astype(complex)
def setup():
 d=thermal_transition_jumps(H,A,beta=BETA,gamma=GAMMA)
 return canonical_gibbs(H,BETA),local_lindblad(P),d['jumps'],d['rates'],d['energies']
def comm(X):
 return 1j*(G@X-X@G)

class ThermalLindbladWard(unittest.TestCase):
 def test_spin_orbit_two_site_and_physical_spin_nonconservation(self):
  self.assertEqual(H.shape,(4,4))
  np.testing.assert_allclose(H,H.conj().T,atol=1e-14)
  spin_z=np.kron(np.eye(2),np.diag([1.,-1.]))
  self.assertGreater(np.max(abs(H@spin_z-spin_z@H)),.1)
  self.assertGreater(np.min(np.diff(np.linalg.eigvalsh(H))),.05)

 def test_complete_onsite_kraus_projectors(self):
  np.testing.assert_allclose(P.sum(axis=0),np.eye(4),atol=1e-14)
  for p in P:
   np.testing.assert_allclose(p@p,p,atol=1e-14)
  self.assertEqual(local_lindblad(P).shape,(2,4,4))

 def test_local_dephasing_preserves_each_site_charge(self):
  rng=np.random.default_rng(29)
  x=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
  _,jumps,*_=setup()
  d=dissipator(x,jumps)
  for p in P:self.assertLess(abs(np.trace(p@d)),1e-13)

 def test_local_dephasing_not_a_finite_temperature_thermalizer(self):
  rho,jumps,*_=setup()
  self.assertGreater(np.max(abs(dissipator(rho,jumps))),.005)
  self.assertGreater(np.trace(H@dissipator(rho,jumps)).real,.01)

 def test_thermal_rates_canonical_detailed_balance(self):
  rho,_,jumps,W,e=setup()
  p=np.exp(-BETA*(e-e.min()));p/=p.sum()
  np.testing.assert_allclose(W*p[None,:],W.T*p[:,None],atol=1e-14)
  self.assertGreater(np.min(W[~np.eye(4,dtype=bool)]),1e-3)
  self.assertEqual(jumps.shape,(12,4,4))

 def test_gibbs_stationary_and_independent_population_master(self):
  rho,_,thermal,W,e=setup()
  self.assertLess(np.max(abs(dissipator(rho,thermal))),1e-13)
  p=np.exp(-BETA*(e-e.min()));p/=p.sum()
  np.testing.assert_allclose(population_rhs(W,p),0,atol=1e-13)
  self.assertAlmostEqual(np.trace(rho).real,1.,places=13)

 def test_kms_bath_global_particle_number_but_not_local_charge(self):
  rho,_,thermal,_,_=setup()
  X=P[0]/2
  source=np.array([np.trace(p@dissipator(X,thermal)) for p in P])
  self.assertGreater(abs(source[0]),.02)
  np.testing.assert_allclose(source.sum(),0,atol=1e-13)

 def test_exact_site_continuity_requires_explicit_thermal_bath_source(self):
  rho,local,thermal,_,_=setup()
  X=P[0]/2
  for jumps in (local,thermal):
   terms=continuity_terms(H,X,P,jumps)
   np.testing.assert_allclose(terms['local_residual'],0,atol=2e-13)
   np.testing.assert_allclose(terms['density_rate']+terms['divergence'],
                              terms['bath_source'],atol=2e-13)
  t=continuity_terms(H,X,P,thermal)
  self.assertGreater(np.max(abs(t['density_rate']+t['divergence'])),.02)

 def test_gksl_trace_and_hermiticity_preservation_both_baths(self):
  rng=np.random.default_rng(77)
  z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
  X=z+z.conj().T
  _,local,thermal,_,_=setup()
  for jumps in (local,thermal):
   d=dissipator(X,jumps)
   self.assertLess(abs(np.trace(d)),2e-13)
   np.testing.assert_allclose(d,d.conj().T,atol=2e-13)

 def test_independent_kronecker_liouvillian_vs_direct_gksl(self):
  rng=np.random.default_rng(55)
  z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
  _,local,thermal,_,_=setup()
  for jumps in (local,thermal):
   direct=master_rhs(H,z,jumps)
   dense=(liouvillian(H,jumps)@z.reshape(16,order='F')).reshape((4,4),order='F')
   np.testing.assert_allclose(dense,direct,atol=2e-13)

 def test_independent_time_evolution_positive_trace_and_rk4(self):
  X=P[0]/2
  _,local,thermal,_,_=setup()
  for jumps in (local,thermal):
   rk=evolved_rk4(H,X,jumps,dt=.01,steps=150)
   exact=evolved_dense(H,X,jumps,time=1.5)
   np.testing.assert_allclose(rk,exact,atol=3e-9)
   self.assertAlmostEqual(np.trace(rk).real,1.,places=11)
   self.assertGreater(np.min(np.linalg.eigvalsh((rk+rk.conj().T)/2)),-2e-10)

 def test_local_u1_joint_transform_covariant_thermal_bath(self):
  rho,_,thermal,_,_=setup()
  X=P[0]/2+rho/3
  phi=.37
  rotated=np.array([gauge_rotate(j,G,phi) for j in thermal])
  rr=gauge_rotate(X,G,phi)
  np.testing.assert_allclose(dissipator(rr,rotated),
                       gauge_rotate(dissipator(X,thermal),G,phi),atol=2e-13)
  np.testing.assert_allclose(master_rhs(gauge_rotate(H,G,phi),rr,rotated),
                       gauge_rotate(master_rhs(H,X,thermal),G,phi),atol=2e-13)

 def test_stale_unrotated_jumps_break_local_gauge_covariance(self):
  rho,_,thermal,_,_=setup()
  X=P[0]/2
  phi=.39
  wrong=dissipator(gauge_rotate(X,G,phi),thermal)
  correct=gauge_rotate(dissipator(X,thermal),G,phi)
  self.assertGreater(np.max(abs(wrong-correct)),.001)

 def test_bath_vertex_is_finite_difference_of_rotated_jump_superoperator(self):
  rho,_,jumps,_,_=setup()
  eps=2e-5
  def change(phi):
   rotated=np.array([gauge_rotate(j,G,phi) for j in jumps])
   return dissipator(rho,rotated)
  finite=(change(eps)-change(-eps))/(2*eps)
  analytic=bath_gauge_vertex(rho,jumps,G)
  np.testing.assert_allclose(finite,analytic,atol=1e-9)

 def test_missing_bath_vertex_spoils_gauge_rotated_gibbs_stationarity(self):
  rho,_,thermal,_,_=setup()
  dr=comm(rho)
  bad=dissipator(dr,thermal)
  proper=bad+bath_gauge_vertex(rho,thermal,G)
  self.assertGreater(np.max(abs(bad)),.005)
  np.testing.assert_allclose(proper,0,atol=2e-13)

 def test_bad_rate_and_gauge_input_fail_closed(self):
  with self.assertRaisesRegex(ValueError,'positive'):canonical_gibbs(H,beta=0)
  with self.assertRaisesRegex(ValueError,'gamma'):local_lindblad(P,gamma=-.1)
  with self.assertRaises(ValueError):thermal_transition_jumps(H,A,beta=-.1)
  with self.assertRaises(ValueError):thermal_transition_jumps(H,np.ones((3,3)))
  with self.assertRaisesRegex(ValueError,'gauge'):gauge_rotate(H,np.ones((4,4)),.2)
  with self.assertRaisesRegex(ValueError,'jumps'):dissipator(H,np.zeros((1,3,3)))
if __name__=='__main__':unittest.main()
