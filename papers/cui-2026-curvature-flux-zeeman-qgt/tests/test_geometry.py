import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from flux_lab.core import (model,qgt_eigen,qgt_projector,qgt_analytic,
  omega_analytic_massive,winding_field,flux_circle,flux_ellipse,winding_ellipse,berry_phase_wilson_circle, flux_circle_from_qgt, divergence_disk_massive)

class Geometry(unittest.TestCase):
 def test_three_routes_for_full_complex_tensor(self):
  for kx,ky,m,vx,vy in [(0.7,-.4,.3,1.,1.),(-.3,.2,0.,1.4,.7),(.21,.53,.41,-.6,1.2)]:
   d,v=model(kx,ky,m,vx,vy)
   for band in [-1,1]:
    a=qgt_analytic(d,v,band)
    for f in [qgt_eigen,qgt_projector]:
     b=f(d,v,band)
     np.testing.assert_allclose(a['T'],b['T'],atol=2e-14)
     for key in ['gN','gA','omegaN','omegaA']:
      np.testing.assert_allclose(a[key],b[key],atol=2e-14)
 def test_gauge_basis_covariance(self):
  d,v=model(.37,-.88,.19)
  theta=.313
  U=np.array([[np.cos(theta),np.sin(theta)*np.exp(1.23j)],[-np.sin(theta)*np.exp(-1.23j),np.cos(theta)]])
  np.testing.assert_allclose(qgt_eigen(d,v,-1)['T'],qgt_eigen(d,v,-1,spin_basis=U)['T'],atol=2e-14)
 def test_hermitian_qgt_limit(self):
  z=np.array([[1.5,0.6-0.3j,0.], [.6+.3j,.8,0.],[0.,0.,0.]])
  self.assertTrue(np.allclose(z,z.conj().T))
  from flux_lab.core import _cross_geom
  terms=_cross_geom(z)
  self.assertLess(np.max(abs(terms['gA'])),1e-14)
  self.assertLess(np.max(abs(terms['omegaA'])),1e-14)
 def test_band_parity(self):
  d,v=model(.8,.3,.2)
  minus=qgt_eigen(d,v,-1);plus=qgt_eigen(d,v,1)
  for k in ['gN','omegaA']: np.testing.assert_allclose(minus[k],plus[k],atol=3e-14)
  for k in ['gA','omegaN']: np.testing.assert_allclose(minus[k],-plus[k],atol=3e-14)
 def test_massive_omega_exact(self):
  for kx,ky,m in [(1.3,.2,.4),(-.2,.4,0.),(.23,-.53,-.13)]:
   d,v=model(kx,ky,m)
   np.testing.assert_allclose(qgt_eigen(d,v,-1)['omegaA'],omega_analytic_massive(kx,ky,m),atol=2e-14)
 def test_planar_hodge_anisotropic(self):
  kx,ky,vx,vy=.42,-.33,1.7,-.8
  d,v=model(kx,ky,0,vx,vy)
  w=winding_field(kx,ky,vx,vy)
  np.testing.assert_allclose(qgt_projector(d,v,1)['omegaA'][:2],[-w[1]/2,w[0]/2],atol=2e-14)
 def test_flux_circle_massless(self):
  for N in [128,512,2048]:
   self.assertLess(abs(flux_circle(1.1,0,N)+np.pi),2e-3)
 def test_mass_regularization_order(self):
  r=.8;m=.16
  for N in [256,1024]:
   self.assertLess(abs(flux_circle(r,m,N)+np.pi*r*r/(r*r+m*m)),3e-4)
  self.assertLess(abs(flux_circle(.005,m,1024)),0.005)
 def test_ellipse_topology_chirality(self):
  for vx,vy,expected in [(1.4,.8,1.),(1.4,-.8,-1.)]:
   w=winding_ellipse(.6,1.8,vx,vy)
   f=flux_ellipse(.6,1.8,vx,vy)
   self.assertLess(abs(w-expected),1e-6)
   self.assertLess(abs(f+np.pi*expected),2e-6)
 def test_wilson_berry_phase(self):
  for mass in [0.,.17]:
   gamma=berry_phase_wilson_circle(1.,mass,points=512,band=-1)
   target=np.pi*(1-mass/np.sqrt(1+mass*mass))
   self.assertLess(abs(np.angle(np.exp(1j*(gamma-target)))),5e-5)
 def test_raw_qgt_contour_vs_analytic_geometry(self):
  for m in [0.,.18]:
   target=-np.pi*.7**2/(.7**2+m*m)
   self.assertLess(abs(flux_circle_from_qgt(.7,m,1024,-1)-target),2e-5)
   self.assertLess(abs(flux_circle_from_qgt(.7,m,1024,1)-target),2e-5)
 def test_finite_mass_gauss_surface(self):
  R=.8;m=.16
  target=-np.pi*R*R/(R*R+m*m)
  self.assertLess(abs(divergence_disk_massive(R,m,10000)-target),2e-7)
  self.assertLess(abs(divergence_disk_massive(R,m,10000)-flux_circle(R,m,4096)),2e-6)
 def test_degenerate_node_guard(self):
  d,v=model(0.,0.,0.)
  with self.assertRaises(ValueError):qgt_eigen(d,v,-1)
  with self.assertRaises(ValueError):qgt_projector(d,v,1)

if __name__=='__main__':unittest.main()
