"""v60 physics gates for matrix principal-layer leads and same-model Chern marker."""
import unittest
from dataclasses import replace
import numpy as np
from scipy.sparse import eye,lil_matrix
from scipy.sparse.linalg import splu
from qgmap.matrix_leads_v60 import (MatrixLead,surface_g_matrix,selfenergy_matrix,matrix_caroli,
 embedded_four_site_leads,integrated_matrix_leads)
from qgmap.qwz_finite_marker_v60 import qwz_open,local_chern_marker,compare_finite_periodic
from qgmap.lead_transport_v57 import diamond_graph
from qgmap.lead_selfenergy_v59 import reference_device,integrate_leads
from qgmap.embedding_magnetic_v55 import peierls_field,field_gauge_unitary

class V60Physics(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.pos,cls.h,_=diamond_graph(1)
  cls.leads=embedded_four_site_leads(cls.h,cls.pos,orbital_matrix=True,eta=.003)
  cls.a=integrated_matrix_leads(cls.h,cls.pos,cls.leads,epsrel=1e-6)

 def test_01_scalar_retarded_analytic(self):
  for E in (-3.,-1.,0.,.4,1.7,3.2):
   t=2.;z=E+.004j;r=np.sqrt(z*z-4*t*t)
   if r.imag<0:r=-r
   g=surface_g_matrix(E,np.array([[0.]]),np.array([[t]]),eta=.004)[0,0]
   self.assertAlmostEqual(g,(z-r)/(2*t*t),delta=1e-9)

 def test_02_matrix_Dyson_noncommuting(self):
  l=self.leads[0];E=.33;eta=.05;g=surface_g_matrix(E,l.h00,l.h01,eta)
  np.testing.assert_allclose(np.linalg.inv(g),(E+1j*eta)*np.eye(2)-l.h00-l.h01@g@l.h01.conj().T,atol=4e-11)
  self.assertGreater(abs(g[0,1]),.02)

 def test_03_independent_160_layer_sparse_resolvent(self):
  l=self.leads[0];N=160;d=2*N;E=.4;eta=.32;mat=lil_matrix((d,d),dtype=complex)
  for j in range(N):
   mat[2*j:2*j+2,2*j:2*j+2]=l.h00
   if j+1<N:
    mat[2*j:2*j+2,2*j+2:2*j+4]=l.h01
    mat[2*j+2:2*j+4,2*j:2*j+2]=l.h01.conj().T
  rhs=np.zeros((d,2),complex);rhs[:2]=np.eye(2)
  finite=splu((E+1j*eta)*eye(d)-mat.tocsc()).solve(rhs)[:2]
  np.testing.assert_allclose(finite,surface_g_matrix(E,l.h00,l.h01,eta),atol=1e-10)

 def test_04_positive_spectral_and_Gamma_matrix(self):
  for l in self.leads:
   sig,_,gamma,_=selfenergy_matrix(.21,l)
   self.assertGreaterEqual(np.linalg.eigvalsh(gamma).min(),-1e-10)
   np.testing.assert_allclose(gamma,1j*(sig-sig.conj().T),atol=1e-13)

 def test_05_orbital_unitary_invariance_in_lead(self):
  l=self.leads[0];A=np.array([[.3+.2j,.4],[.1j,.2-.5j]])
  U,_=np.linalg.qr(A)
  b=MatrixLead(U.conj().T@l.h00@U,U.conj().T@l.h01@U,l.coupling@U,l.eta)
  for E in (-.4,.13,1.7):
   x=selfenergy_matrix(E,l);y=selfenergy_matrix(E,b)
   for a,c in zip(x,y):np.testing.assert_allclose(a,c,atol=1e-10)

 def test_06_nonlocal_lead_contact_covariant_Sigma_derivative(self):
  l=self.leads[0];c=l.coupling.copy();c[1]=[.11,.22j];l=replace(l,coupling=c)
  G,dG=field_gauge_unitary(self.pos+np.array([.3,.5]),.31,-.87)
  x=selfenergy_matrix(.24,l);y=selfenergy_matrix(.24,replace(l,coupling=G@c),dG@c)
  np.testing.assert_allclose(y[0],G@x[0]@G.conj().T,atol=1e-12)
  np.testing.assert_allclose(y[1],dG@x[0]@G.conj().T+G@x[0]@dG.conj().T,atol=1e-12)
  self.assertGreater(abs(y[2][0,1]),.001)

 def test_07_nonlocal_Caroli_full_matrix_trace(self):
  p=peierls_field(self.h,self.pos,.31,-.87)
  ls=[]
  for j,l in enumerate(self.leads):
   c=l.coupling.copy();c[(j+1)%4]=[.06,.09j];ls.append(replace(l,coupling=c))
  T,dT=matrix_caroli(.18,p.hamiltonian,p.magnetic_vertex,ls)
  sigs=[selfenergy_matrix(.18,l) for l in ls]
  G=np.linalg.inv(.18*np.eye(4)-p.hamiltonian-sum((x[0] for x in sigs),np.zeros_like(self.h)))
  for a in range(4):
   for b in range(4):
    if a!=b:self.assertAlmostEqual(T[a,b],np.trace(sigs[a][2]@G@sigs[b][2]@G.conj().T).real,delta=1e-12)
  self.assertGreater(np.max(abs(dT)),.03)

 def test_08_full_magnetic_derivative_at_fixed_energy(self):
  p=peierls_field(self.h,self.pos,.31,-.87)
  _,dT=matrix_caroli(.18,p.hamiltonian,p.magnetic_vertex,self.leads)
  def f(B):
   x=peierls_field(self.h,self.pos,B,-.87)
   return matrix_caroli(.18,x.hamiltonian,x.magnetic_vertex,self.leads)[0]
  np.testing.assert_allclose(dT,(f(.31001)-f(.30999))/.00002,atol=3e-9)

 def test_09_contact_B_gauge_derivative_is_essential(self):
  pos=self.pos+np.array([.3,.5]);G,dG=field_gauge_unitary(pos,.31,-.87)
  ps=peierls_field(self.h,pos,.31,-.87)
  pl=peierls_field(self.h,pos,.31,-.87,gauge='landau_y')
  ls=[]
  for j,l in enumerate(self.leads):
   c=l.coupling.copy();c[(j+1)%4]=[.06,.09j];ls.append(replace(l,coupling=c))
  L=[replace(l,coupling=G@l.coupling) for l in ls];dl=[dG@l.coupling for l in ls]
  a,b=matrix_caroli(.25,ps.hamiltonian,ps.magnetic_vertex,ls)
  c,d=matrix_caroli(.25,pl.hamiltonian,pl.magnetic_vertex,L,dl)
  np.testing.assert_allclose(a,c,atol=3e-11);np.testing.assert_allclose(b,d,atol=4e-10)
  self.assertGreater(np.max(abs(matrix_caroli(.25,pl.hamiltonian,pl.magnetic_vertex,L)[1]-b)),1e-4)

 def test_10_magnetic_gauge_for_whole_fourprobe(self):
  pos=self.pos+np.array([.3,.5])
  a=integrated_matrix_leads(self.h,pos,self.leads,epsrel=1e-6)
  b=integrated_matrix_leads(self.h,pos,self.leads,epsrel=1e-6,gauge='landau_y')
  np.testing.assert_allclose(a['T'],b['T'],atol=8e-9)
  np.testing.assert_allclose(a['dT_dB'],b['dT_dB'],atol=8e-8)
  self.assertAlmostEqual(a['observables']['dHall_dB'],b['observables']['dHall_dB'],delta=2e-8)

 def test_11_complete_finiteB_recomputed_probe(self):
  def f(B):return integrated_matrix_leads(self.h,self.pos,self.leads,field=B,epsrel=1e-6)['observables']['hall_ratio']
  self.assertAlmostEqual((f(.3101)-f(.3099))/.0002,self.a['observables']['dHall_dB'],delta=2e-6)

 def test_12_scalar_eta_small_recovers_v59(self):
  l=embedded_four_site_leads(self.h,self.pos,orbital_matrix=False,eta=.001)
  a=integrated_matrix_leads(self.h,self.pos,l,epsrel=1e-6)
  h,p,ls=reference_device();b=integrate_leads(h,p,ls)
  self.assertAlmostEqual(a['observables']['hall_ratio'],b['observables']['hall_ratio'],delta=1e-5)
  self.assertAlmostEqual(a['observables']['dHall_dB'],b['observables']['dHall_dB'],delta=2e-5)

 def test_13_eta_convergence_not_silently_set_to_zero(self):
  ys=[]
  for eta in (.01,.003,.001):
   l=embedded_four_site_leads(self.h,self.pos,orbital_matrix=False,eta=eta)
   ys.append(integrated_matrix_leads(self.h,self.pos,l,epsrel=2e-6)['observables']['hall_ratio'])
  self.assertLess(abs(ys[2]-ys[1]),abs(ys[1]-ys[0]))
  self.assertGreater(abs(ys[2]-ys[1]),1e-8)

 def test_14_matrix_lead_bad_inputs_fail_closed(self):
  l=self.leads[0]
  for call in (lambda:MatrixLead(l.h00,l.h01,np.ones((4,1))),
               lambda:MatrixLead(l.h00,l.h01,l.coupling,eta=0),
               lambda:surface_g_matrix(.2,l.h00,l.h01,eta=-.1),
               lambda:selfenergy_matrix(.1,l,np.ones((2,2))),
               lambda:matrix_caroli(.3,self.h,np.eye(4),self.leads[:1])):
   with self.assertRaises(ValueError):call()

 def test_15_open_lattice_Fourier_matches_same_periodic_QWZ(self):
  from qgmap.periodic_kubo_v59 import qwz_h,sz
  h,_=qwz_open(5,5,-1)
  i=2*(1*5+1);Tx=h[i:i+2,i+2:i+4];Ty=h[i:i+2,i+10:i+12]
  for x,y in ((.27,.36),(-1.,2.)):
   v=-sz+Tx*np.exp(1j*x)+Tx.conj().T*np.exp(-1j*x)+Ty*np.exp(1j*y)+Ty.conj().T*np.exp(-1j*y)
   np.testing.assert_allclose(v,qwz_h(x,y,-1),atol=2e-15)

 def test_16_finite_marker_trace_vanishes_edge_compensates(self):
  h,xy=qwz_open(9,9,-1);m=local_chern_marker(h,xy)
  self.assertLess(abs(m['sum_over_all_orbitals']),1e-10)
  self.assertGreater(m['marker_site'][4*9+4],.99)
  self.assertLess(m['marker_site'][0],0)

 def test_17_same_model_center_converges_to_periodic_chern(self):
  vals=[compare_finite_periodic(L,-1,max(1,L//4)) for L in (5,7,9,11)]
  self.assertLess(abs(vals[-1]['center_marker']-1),abs(vals[0]['center_marker']-1))
  self.assertLess(abs(vals[-1]['interior_marker_mean']-1),.004)
  self.assertAlmostEqual(vals[-1]['periodic_chern'],1,delta=1e-12)

 def test_18_correct_occupied_chern_sign_and_trivial_phase(self):
  for mass,target in ((-3,0),(-1,1),(1,-1),(3,0)):
   x=compare_finite_periodic(9,mass,2)
   self.assertAlmostEqual(x['center_marker'],target,delta=.045)
   self.assertAlmostEqual(x['periodic_chern'],target,delta=1e-12)

 def test_19_conjugation_reverses_finite_chern_sign(self):
  h,xy=qwz_open(7,7,-1)
  a=local_chern_marker(h,xy)['marker_site']
  b=local_chern_marker(h.conj(),xy)['marker_site']
  np.testing.assert_allclose(a,-b,atol=2e-12)

 def test_20_orbital_basis_rotation_preserves_site_sum(self):
  h,xy=qwz_open(5,5,-1);A=np.array([[.3+.2j,.4],[.1j,.2-.5j]]);U,_=np.linalg.qr(A)
  W=np.kron(np.eye(25),U)
  a=local_chern_marker(h,xy)['marker_site']
  b=local_chern_marker(W.conj().T@h@W,xy)['marker_site']
  np.testing.assert_allclose(a,b,atol=3e-12)

 def test_21_occupied_geometry_not_certified_transport(self):
  a=compare_finite_periodic(7,-1,1)
  self.assertIn('NOT_TRANSPORT',a['status'])
  self.assertNotIn('conductance',a)
  with self.assertRaises(ValueError):qwz_open(2,3,-1)
  with self.assertRaises(ValueError):compare_finite_periodic(7,margin=4)

if __name__=='__main__':unittest.main()
