"""Round30 exact fermionic CAR, interacting grand Gibbs, KMS and source Ward."""
import unittest
import numpy as np
from flux_lab.fermionic_fock_bath import (
 fermion_ops,lift_one_body,fock_soc_hubbard,grand_gibbs,
 thermal_sector_jumps,local_charge_balance,unitary_phase,bath_contact,
 sector_gibbs_weights)
from flux_lab.thermal_lindblad_ward import dissipator,liouvillian,master_rhs

BETA=.83;MU=.16;GAMMA=.14

def setup():
 d=fock_soc_hubbard()
 t=thermal_sector_jumps(d['H'],d['N'],d['coupling'],beta=BETA,gamma=GAMMA)
 rho=grand_gibbs(d['H'],d['N'],beta=BETA,mu=MU)
 return d,t,rho

class InteractingFockThermalBath(unittest.TestCase):
 def test_four_mode_car_for_all_pairs(self):
  c,cd=fermion_ops()
  for i in range(4):
   for j in range(4):
    np.testing.assert_allclose(c[i]@cd[j]+cd[j]@c[i],np.eye(16) if i==j else 0,atol=1e-14)
    np.testing.assert_allclose(c[i]@c[j]+c[j]@c[i],0,atol=1e-14)
    np.testing.assert_allclose(cd[i]@cd[j]+cd[j]@cd[i],0,atol=1e-14)

 def test_pauli_blocked_fermion_hop_and_occupied_target(self):
  c,cd=fermion_ops();hop=cd[1]@c[0]
  blocked=np.eye(16)[:,0b0011]
  available=np.eye(16)[:,0b0001]
  np.testing.assert_allclose(hop@blocked,0,atol=1e-14)
  np.testing.assert_allclose(hop@available,np.eye(16)[:,0b0010],atol=1e-14)
  self.assertAlmostEqual(np.linalg.norm(hop@available),1)

 def test_interacting_hubbard_four_mode_charge_and_physical_spin(self):
  d,_,_=setup()
  self.assertEqual(d['H'].shape,(16,16))
  np.testing.assert_allclose(d['H'],d['H'].conj().T,atol=1e-14)
  np.testing.assert_allclose(d['H']@d['N']-d['N']@d['H'],0,atol=2e-14)
  self.assertGreater(np.max(abs(d['H']@d['zeeman']-d['zeeman']@d['H'])),.1)
  np.testing.assert_allclose(d['occupations'].sum(axis=0),d['N'],atol=1e-14)

 def test_exact_sector_dimensions_and_number_conserving_jumps(self):
  d,t,_=setup()
  self.assertEqual([int(np.sum(t['numbers']==k)) for k in range(5)],[1,4,6,4,1])
  for L in t['jumps']:
   self.assertLess(np.max(abs(L@d['N']-d['N']@L)),2e-13)
  self.assertGreater(len(t['jumps']),18)

 def test_fermionic_grand_gibbs_and_occupation_range(self):
  d,t,rho=setup()
  self.assertAlmostEqual(np.trace(rho).real,1.,places=13)
  self.assertGreater(np.min(np.linalg.eigvalsh(rho)),0)
  fills=[np.trace(n@rho).real for n in d['occupations']]
  self.assertTrue(all(0<x<1 for x in fills))
  self.assertTrue(0<np.trace(d['N']@rho).real<4)

 def test_noninteracting_fock_occ_matches_independent_fermi_formula(self):
  d=fock_soc_hubbard(interaction=0,cross=0)
  rho=grand_gibbs(d['H'],d['N'],beta=BETA,mu=MU)
  e,U=np.linalg.eigh(d['h'])
  f=1/(1+np.exp(BETA*(e-MU)))
  orbital_density=(U*f)@U.conj().T
  exact=np.array([np.trace(ni@rho).real for ni in d['occupations']])
  np.testing.assert_allclose(exact,np.diag(orbital_density).real,atol=1e-12)

 def test_interactions_measurably_change_orbital_occupation(self):
  d,t,rho=setup()
  reference=fock_soc_hubbard(interaction=0,cross=0)
  bare=grand_gibbs(reference['H'],reference['N'],beta=BETA,mu=MU)
  occupied=np.array([np.trace(n@rho).real for n in d['occupations']])
  uncorrelated=np.array([np.trace(n@bare).real for n in reference['occupations']])
  self.assertGreater(np.max(abs(occupied-uncorrelated)),.01)

 def test_grand_canonical_kms_pairwise_balanced_with_sector_constraint(self):
  _,t,_=setup()
  E,N,W=t['energies'],t['numbers'],t['rates']
  p=sector_gibbs_weights(E,N,beta=BETA,mu=MU)
  np.testing.assert_allclose(W*p[None,:],W.T*p[:,None],atol=2e-14)
  self.assertEqual(np.count_nonzero(W[N[:,None]!=N[None,:]]),0)
  self.assertGreater(np.max(W),.01)

 def test_exact_interacting_grand_gibbs_stationary_for_chosen_thermal_bath(self):
  d,t,rho=setup()
  np.testing.assert_allclose(dissipator(rho,t['jumps']),0,atol=1e-12)
  p=sector_gibbs_weights(t['energies'],t['numbers'],beta=BETA,mu=MU)
  independent=(t['eigenvectors']*p)@t['eigenvectors'].conj().T
  np.testing.assert_allclose(rho,independent,atol=2e-12)

 def test_number_resolved_bath_cannot_fix_chemical_potential(self):
  d,t,_=setup()
  alternate=grand_gibbs(d['H'],d['N'],beta=BETA,mu=-.62)
  np.testing.assert_allclose(dissipator(alternate,t['jumps']),0,atol=1e-12)
  original=grand_gibbs(d['H'],d['N'],beta=BETA,mu=MU)
  self.assertGreater(abs(np.trace(d['N']@(original-alternate))),.1)

 def test_number_conservation_for_arbitrary_fock_density(self):
  d,t,_=setup()
  rng=np.random.default_rng(30);z=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
  X=z@z.conj().T;X/=np.trace(X)
  diss=dissipator(X,t['jumps'])
  self.assertLess(abs(np.trace(diss)),2e-13)
  self.assertLess(abs(np.trace(d['N']@diss)),2e-13)
  self.assertLess(abs(np.trace(d['N']@d['N']@diss)),2e-13)

 def test_local_site_ward_has_nonzero_bath_source(self):
  d,t,_=setup()
  X=np.zeros((16,16),complex);X[0b0011,0b0011]=1.
  a=local_charge_balance(d,X,t['jumps'])
  np.testing.assert_allclose(a['continuity'],0,atol=2e-13)
  self.assertLess(abs(a['total_source']),2e-13)
  self.assertGreater(np.max(abs(a['bath_source'])),.002)
  self.assertGreater(np.max(abs(a['density_rate']+np.array([a['current'],-a['current']]))),.002)

 def test_independent_kronecker_liouville_matches_fock_jump_master(self):
  d,t,_=setup()
  X=np.zeros((16,16),complex);X[0b0011,0b0011]=1.
  independent=(liouvillian(d['H'],t['jumps'])@X.reshape(-1,order='F')).reshape((16,16),order='F')
  np.testing.assert_allclose(independent,master_rhs(d['H'],X,t['jumps']),atol=4e-13)

 def test_local_manybody_u1_jump_covariance_and_stale_control(self):
  d,t,rho=setup();G=d['gauge'];phi=.41
  X=np.zeros((16,16),complex);X[0b0011,0b0011]=.6;X+=.4*rho
  rot=lambda x:unitary_phase(x,G,phi)
  jumps=np.array([rot(j) for j in t['jumps']])
  np.testing.assert_allclose(dissipator(rot(X),jumps),rot(dissipator(X,t['jumps'])),atol=3e-13)
  wrong=dissipator(rot(X),t['jumps'])
  self.assertGreater(np.max(abs(wrong-rot(dissipator(X,t['jumps'])))),.001)

 def test_bath_gauge_vertex_independent_finite_amplitude_derivative(self):
  d,t,rho=setup();G=d['gauge']
  eps=2e-5
  def Dphi(phi):
   js=np.array([unitary_phase(j,G,phi) for j in t['jumps']])
   return dissipator(rho,js)
  numeric=(Dphi(eps)-Dphi(-eps))/(2*eps)
  analytic=bath_contact(rho,t['jumps'],G)
  np.testing.assert_allclose(numeric,analytic,atol=2e-9)
  dr=1j*(G@rho-rho@G)
  np.testing.assert_allclose(dissipator(dr,t['jumps'])+analytic,0,atol=2e-12)
  self.assertGreater(np.max(abs(dissipator(dr,t['jumps']))),1e-4)

 def test_magnetic_bath_vertex_is_required_for_thermal_stationarity(self):
  d,t,rho=setup()
  eps=3e-5;H=d['H'];V=d['zeeman']
  def at(sign):
   hp=H+sign*eps*V
   rhop=grand_gibbs(hp,d['N'],beta=BETA,mu=MU)
   j=thermal_sector_jumps(hp,d['N'],d['coupling'],beta=BETA,gamma=GAMMA)['jumps']
   return rhop,j
  rp,jp=at(1);rm,jm=at(-1)
  dr=(rp-rm)/(2*eps)
  vertex=(dissipator(rho,jp)-dissipator(rho,jm))/(2*eps)
  state_term=dissipator(dr,t['jumps'])
  self.assertGreater(np.max(abs(state_term)),1e-4)
  np.testing.assert_allclose(state_term+vertex,0,atol=3e-8)
  self.assertGreater(np.max(abs(rp-rm))/(2*eps),.01)

 def test_invalid_thermal_fermionic_data_fail_closed(self):
  d,t,_=setup()
  with self.assertRaises(ValueError):fermion_ops(3)
  with self.assertRaises(ValueError):grand_gibbs(d['H'],d['N'],beta=0)
  with self.assertRaises(ValueError):thermal_sector_jumps(d['H'],d['N'],d['coupling'],gamma=-1)
  with self.assertRaisesRegex(ValueError,'conserve total N'):
   thermal_sector_jumps(d['H']+d['annihilate'][0]+d['create'][0],d['N'],d['coupling'])
  with self.assertRaisesRegex(ValueError,'charge gauge'):
   unitary_phase(d['H'],np.ones((16,16)),.2)

if __name__=='__main__':unittest.main()
