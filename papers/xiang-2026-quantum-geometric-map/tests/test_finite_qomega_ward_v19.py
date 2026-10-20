"""Original v19 finite-qomega resolvent Ward method tests, not full PRB acceptance."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.finite_qomega_ward_v19 import finite_ward,two_energy_greens
from qgmap.warped_disorder_v17 import Lattice,k_grid,h_and_vel
from qgmap.microscopic_ward_v11 import I,ladder_vertex,scba,greens

class TestFiniteQOmegaWardV19(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=Lattice(n=41,alpha=.16,lam=.6,eta=.08)

    def good(self,**kwargs):
        r=finite_ward(self.p,**kwargs)
        self.assertLess(r['residuals']['wti'],1e-9)
        self.assertLess(r['residuals']['resolvent'],1e-12)
        self.assertLess(r['residuals']['ladder'],1e-10)
        self.assertLess(r['residuals']['selfenergy'],1e-10)
        return r

    def test_finite_omega_qy(self):self.good(q=.08,omega=.04)
    def test_finite_omega_qx(self):self.good(q=.11,omega=.03,direction=0)
    def test_zero_q_finite_omega(self):self.good(q=0.,omega=.03)
    def test_finite_q_static(self):self.good(q=.08,omega=0.)
    def test_both_zero_finite_eta(self):self.good(q=0.,omega=0.)
    def test_shift_protocol(self):self.good(q=.08,omega=.04,b=.025,protocol='shift')
    def test_onsite_magnetic_field(self):self.good(q=.08,omega=.04,b=.025,protocol='onsite')

    def test_disorder_free_resolvent(self):
        r=finite_ward(replace(self.p,disorder_variance=0.),q=.12,omega=.02)
        self.assertLess(r['residuals']['wti'],2e-13)
        self.assertLess(r['residuals']['bare'],1e-14)
        self.assertLess(r['density_ladder_correction_norm'],1e-14)

    def test_bare_vertex_fails_nonzero_disorder(self):
        r=self.good(q=.08,omega=.04)
        self.assertGreater(r['residuals']['bare'],1e-3)
        self.assertGreater(r['residuals']['bare']/r['residuals']['wti'],1e8)

    def test_external_eta_not_physical_continuity(self):
        r=self.good(q=.08,omega=.04)
        self.assertAlmostEqual(r['residuals']['physical_frequency_without_regulator'],
           r['residuals']['missing_regulator_expected'],delta=2e-10)
        self.assertGreater(r['residuals']['physical_frequency_without_regulator'],.05)

    def test_density_vertex_q0_matches_v17(self):
        from qgmap.warped_disorder_v17 import density_ward
        a=density_ward(self.p)
        b=finite_ward(self.p,q=0.,omega=0.)
        self.assertLess(a['abs_error'],2e-10)
        self.assertLess(b['residuals']['wti'],2e-10)

    def test_nonzero_q_bare_midpoint_current_is_inexact(self):
        q=.7
        k=k_grid(self.p)
        upper=k.copy();upper[...,1]+=q/2
        lower=k.copy();lower[...,1]-=q/2
        hp,_,_=h_and_vel(upper,self.p)
        hm,_,_=h_and_vel(lower,self.p)
        _,_,v=h_and_vel(k,self.p)
        err=np.max(abs(hp-hm-q*v))
        self.assertGreater(err,1e-3)
        self.assertLess(err,.5)

    def test_density_near_singular_channel_sensitive_to_order(self):
        eta=.005
        p=replace(self.p,eta=eta,n=61)
        zero=finite_ward(p,q=0.,omega=0.)
        freq=finite_ward(p,q=0.,omega=.025)
        momentum=finite_ward(p,q=.06,omega=0.)
        self.assertGreater(zero['density_vertex_norm'],2*freq['density_vertex_norm'])
        self.assertGreater(zero['density_vertex_norm'],2*momentum['density_vertex_norm'])
        self.assertLess(max(zero['residuals']['wti'],freq['residuals']['wti'],momentum['residuals']['wti']),1e-9)

    def test_constant_basis_unitary_covariance_of_ladders(self):
        p=self.p
        x=two_energy_greens(p,qy=.085,omega=.039)
        q=.085;z=x['regulator'];phi=.37
        sx=np.array([[0,1],[1,0]],complex)
        U=np.cos(phi)*I+1j*np.sin(phi)*sx
        ud=U.conj().T
        rotate=lambda x:ud@x@U
        d=(x['hp']-x['hm'])/q
        for v in (np.broadcast_to(I,d.shape),d):
            v0,c0,r0=ladder_vertex(x['gr'],x['ga'],v,p.disorder_variance)
            vr,cr,rr=ladder_vertex(rotate(x['gr']),rotate(x['ga']),rotate(v),p.disorder_variance)
            self.assertLess(np.max(abs(rotate(v0)-vr)),2e-12)
            self.assertLess(np.max(abs(rotate(c0)-cr)),2e-12)

    def test_bad_q_direction_rejected(self):
        with self.assertRaises(ValueError):finite_ward(self.p,q=.05,direction=4)
        with self.assertRaises(ValueError):finite_ward(self.p,q=np.inf)

if __name__=='__main__':unittest.main()
