"""Fermi-sea Hall consistency gates separate from symmetric Kubo and driven EB."""
import unittest
import numpy as np
from qgmap.streda_clean_sea_v38 import (
    qwz_model,bulk_sea,sea_per_k,direct_energy_sea,clean_berry,
    fhs_chern,sea_matrix_integrand,symmetric_spectral_xx,metallic_pole_audit)

class CleanFermiSea(unittest.TestCase):
    def test_qwz_independent_topology(self):
        for m,cn in [(-3,0),(-1,-1),(1,1),(3,0)]:
            h,vx,vy=qwz_model(45,m)
            self.assertAlmostEqual(fhs_chern(h),cn,delta=1e-10)
            bare=bulk_sea(45,m,eta=.012)
            self.assertLess(abs(bare['occupied_berry_chern']-cn),4e-7)
    def test_green_sea_recovers_integer_in_clean_limit(self):
        for m,cn in [(-1,-1),(1,1)]:
            vals=[bulk_sea(57,m,eta=e)['sea_sigma_over_e2h'] for e in (.16,.06,.02)]
            self.assertLess(abs(vals[-1]-cn),abs(vals[0]-cn))
            self.assertLess(abs(vals[-1]-cn),.03)
    def test_energy_quadrature_independent_of_primitive(self):
        h,vx,vy=qwz_model(23,-1)
        for i,j in ((3,5),(10,12),(17,21)):
            got,err=direct_energy_sea(h[i,j],vx[i,j],vy[i,j],mu=0,eta=.14,lower_cut=-100)
            expected=sea_per_k(h[i,j],vx[i,j],vy[i,j],mu=0,eta=.14)
            self.assertLess(abs(got-expected),1e-7)
            self.assertLess(err,1e-7)
    def test_antisymmetric_xy_bare_sea(self):
        h,vx,vy=qwz_model(21,-1)
        np.testing.assert_allclose(sea_per_k(h,vx,vy),-sea_per_k(h,vy,vx),atol=5e-13)
        self.assertAlmostEqual(sea_matrix_integrand(.23,h[2,7],vx[2,7],vy[2,7],.1),
            -sea_matrix_integrand(.23,h[2,7],vy[2,7],vx[2,7],.1),delta=1e-10)
    def test_u1_gauge_and_constant_su2_basis(self):
        h,vx,vy=qwz_model(41,-1)
        rng=np.random.default_rng(51)
        p=rng.uniform(-np.pi,np.pi,(41,41))
        self.assertAlmostEqual(fhs_chern(h,phases=p),fhs_chern(h),delta=1e-11)
        theta=.43; sx=np.array([[0,1],[1,0]],complex);U=np.cos(theta)*np.eye(2)+1j*np.sin(theta)*sx
        rotated=U.conj().T@h@U
        self.assertAlmostEqual(fhs_chern(rotated),fhs_chern(h),delta=1e-11)
        self.assertLess(np.max(abs(sea_per_k(rotated,U.conj().T@vx@U,U.conj().T@vy@U)-sea_per_k(h,vx,vy))),1e-11)
    def test_metal_finite_eta_not_integer(self):
        r=bulk_sea(51,-1,mu=1.1,eta=.02)
        self.assertGreater(abs(r['occupied_berry_chern']+1),.1)
        # The raw sea integral alone contains a Fermi-level 1/z pole and
        # differs strongly from the occupied-band Berry value in a metal.
        # Removing that pole is ONLY an algebraic diagnostic, not full Streda I.
        self.assertGreater(abs(r['sea_sigma_over_e2h']-r['occupied_berry_chern']),.5)
        d=metallic_pole_audit(101,-1,mu=1.1,eta=.02)
        self.assertLess(abs(d['sea_plus_boundary_residual']),1e-12)
        self.assertLess(abs(d['log_occupancy']-d['occupied_berry_clean']),.06)
    def test_symmetric_spectral_control_vanishes_without_sea(self):
        h,vx,vy=qwz_model(43,-1)
        value1=symmetric_spectral_xx(h,vx,eta=.12)
        value2=symmetric_spectral_xx(h,vx,eta=.04)
        self.assertGreater(value1,0)
        self.assertLess(value2,value1*.2)
        self.assertGreater(abs(bulk_sea(43,-1,eta=.04)['sea_sigma_over_e2h']),.9)
    def test_fail_closed_degenerate(self):
        for mass in (-2,0,2):
            with self.assertRaises(ValueError): qwz_model(mass=mass)
        h,vx,vy=qwz_model(21,-1)
        with self.assertRaises(ValueError):sea_per_k(h,vx,vy,eta=0)
        with self.assertRaises(ValueError):direct_energy_sea(h[0,0],vx[0,0],vy[0,0],lower_cut=2)
        with self.assertRaises(ValueError):fhs_chern(h,phases=np.ones((4,4)))
if __name__=='__main__':unittest.main()
