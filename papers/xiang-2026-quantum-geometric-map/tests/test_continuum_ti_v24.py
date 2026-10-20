"""v24 original tilted Dirac continuum: model-only, no physical DC acceptance."""
from __future__ import annotations
from dataclasses import replace
import unittest
import numpy as np
from qgmap.continuum_ti_v24 import (DiracSpec,quadrature,continuum_matrices,
    analytic_projectors,projector_channels,eigen_channels,eq20_direct,angular_exact)

class TestContinuumV24(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=DiracSpec(nr=112,nphi=160,eta_e=.023,eta_b=.017,omega_e=.016,omega_b=.007)
        cls.project=projector_channels(cls.p)
        cls.eigen=eigen_channels(cls.p)
    def test_01_quadrature_area(self):
        k,w=quadrature(self.p.nr,self.p.nphi,self.p.cutoff)
        self.assertAlmostEqual(np.sum(w),self.p.cutoff**2/(4*np.pi),places=12)
    def test_02_main_h_and_hermiticity(self):
        k,w=quadrature(self.p.nr,self.p.nphi,self.p.cutoff)
        H,Vx,Vy,Z=continuum_matrices(self.p,k)
        self.assertLess(np.max(abs(H-H.conj().transpose(0,2,1))),1e-13)
        self.assertLess(np.max(abs(Vx-Vx.conj().T)),1e-13)
        self.assertLess(np.max(abs(Vy-Vy.conj().T)),1e-13)
    def test_03_two_algorithms_B_first(self):
        self.assertLess(abs(self.project['b_then_e']-self.eigen['b_then_e']),1e-10)
    def test_04_two_algorithms_E_first(self):
        self.assertLess(abs(self.project['e_then_b']-self.eigen['e_then_b']),1e-10)
    def test_05_both_orderings_contribute(self):
        self.assertLess(abs(self.project['total']-self.project['b_then_e']-self.project['e_then_b']),1e-12)
        self.assertGreater(abs(self.project['b_then_e']),1e-5)
        self.assertGreater(abs(self.project['e_then_b']),1e-5)
    def test_06_analytic_projector_derivative(self):
        k,_=quadrature(12,20,1.8)
        step=1e-5
        _,_,rho,dr,_=analytic_projectors(self.p,k)
        e1,p1,r1,_,_=analytic_projectors(self.p,k+np.array([0,step]))
        e2,p2,r2,_,_=analytic_projectors(self.p,k-np.array([0,step]))
        self.assertLess(np.max(abs((r1-r2)/(2*step)-dr)),3e-8)
    def test_07_eq20_angular_reference(self):
        q=replace(self.p,nr=144,nphi=176)
        self.assertLess(abs(eq20_direct(q)-angular_exact(q)),2e-11)
    def test_08_helicity_changes_sign(self):
        q=replace(self.p,omega_e=0,omega_b=0)
        x=eq20_direct(q); y=eq20_direct(replace(q,helicity='sm'))
        self.assertLess(abs(x+y),1e-12)
    def test_09_chemical_potential_changes_sign(self):
        q=replace(self.p,omega_e=0,omega_b=0)
        self.assertLess(abs(eq20_direct(q)+eq20_direct(replace(q,mu=-q.mu))),1e-12)
    def test_10_tilt_changes_sign(self):
        q=replace(self.p,omega_e=0,omega_b=0)
        self.assertLess(abs(eq20_direct(q)+eq20_direct(replace(q,beta=-q.beta))),1e-12)
    def test_11_zero_tilt_and_mu_vanish(self):
        q=replace(self.p,omega_e=0,omega_b=0)
        self.assertLess(abs(eq20_direct(replace(q,beta=0.))),1e-11)
        self.assertLess(abs(eq20_direct(replace(q,mu=0.))),1e-11)
    def test_12_regulator_convergence(self):
        q=DiracSpec(nr=120,nphi=160)
        pole=eq20_direct(q)
        big=.02*projector_channels(replace(q,eta_e=.02,eta_b=.02))['total']
        small=.005*projector_channels(replace(q,eta_e=.005,eta_b=.005))['total']
        self.assertLess(abs(small-pole),9e-6)
        self.assertLess(abs(small-pole),abs(big-pole)*.11)
    def test_13_fixed_frequency_no_pole(self):
        q=replace(self.p,omega_e=.22,omega_b=.07)
        big=.04*projector_channels(replace(q,eta_e=.04,eta_b=.04))['total']
        small=.01*projector_channels(replace(q,eta_e=.01,eta_b=.01))['total']
        self.assertLess(abs(small),abs(big)*.45)
    def test_14_uv_stability(self):
        q=replace(self.p,nr=144,omega_e=0,omega_b=0)
        x=eq20_direct(replace(q,cutoff=2.5));y=eq20_direct(replace(q,cutoff=4.))
        self.assertLess(abs(x-y),1e-10)
    def test_15_T0_closure_vs_main_eq11(self):
        # Eq. 11 reports per physical B, while Eq. 20 couples bbar=-g*mu_B*B.
        q=replace(self.p,nr=144,nphi=176)
        T0=angular_exact(q)/np.tanh(q.mu/(2*q.temp))
        p_eq11=(q.beta**2+2*np.sqrt(1-q.beta**2)-2)/(4*np.pi*q.beta**3)
        self.assertLess(abs(-T0-p_eq11),1e-12)
    def test_16_fail_closed(self):
        for kwargs in ({'beta':1.},{'temp':0},{'eta_e':0},{'helicity':'other'},{'nr':1}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):DiracSpec(**kwargs)

if __name__=='__main__':unittest.main()
