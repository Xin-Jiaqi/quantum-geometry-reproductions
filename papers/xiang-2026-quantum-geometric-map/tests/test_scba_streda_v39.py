"""Model gate: SCBA frequency derivative, RR Ward and exact finite-disorder checks.

These tests are method-level; none proves full disordered Kubo-Streda / EB PHE.
"""
import unittest
import numpy as np
from qgmap.scba_streda_v39 import (
    SCBAConfig, solve_self_energy, derivative_vertex, derivative_audit,
    sea_energy_scan, exact_finite_disorder_sea, at_energy)
from qgmap.streda_clean_sea_v38 import qwz_model, sea_per_k, PAULI

class SCBASeaV39(unittest.TestCase):
    def test_scba_density_ward_RR(self):
        H,vx,vy=qwz_model(21,-1)
        cfg=SCBAConfig(mesh=21,disorder_variance=.12,eta=.12)
        s=at_energy(H,vx,vy,.7,cfg)
        self.assertLess(s['scba_residual'],1e-10)
        self.assertLess(s['vertex_residual'],1e-11)
        self.assertGreater(s['gamma_norm'],.2)
        # Retarded self-energy must have nonpositive anti-Hermitian part.
        sig=s['sigma']
        self.assertLessEqual(float(np.max(np.linalg.eigvalsh((sig-sig.conj().T)/(2j)))),1e-10)

    def test_zero_disorder_RR_vertex_is_identity(self):
        d=derivative_audit(mesh=17,W=0.,h=3e-5)
        self.assertEqual(d['gamma_minus_identity'],0.)
        self.assertEqual(d['derivative_max_abs'],0.)
        self.assertAlmostEqual(d['sample_sea_kernel'],d['sample_sea_naive'],delta=1e-14)

    def test_finite_difference_selfenergy_at_two_disorder_strengths(self):
        for w in (.06,.12):
            with self.subTest(W=w):
                d=derivative_audit(mesh=21,W=w,h=3e-5)
                self.assertLess(d['derivative_max_abs'],5e-7)
                self.assertLess(d['propagator_derivative_max_abs'],5e-6)
                self.assertLess(d['derivative_vertex_residual'],1e-11)

    def test_RR_ward_is_not_constant_sigma_approximation(self):
        d=derivative_audit(mesh=21,W=.12,energy=.7)
        self.assertGreater(d['gamma_minus_identity'],.30)
        self.assertGreater(abs(d['sample_sea_kernel']-d['sample_sea_naive']),.10)

    def test_resolvent_derivative_at_one_energy_by_finite_difference(self):
        H,vx,vy=qwz_model(17,-1)
        cfg=SCBAConfig(mesh=17,disorder_variance=.06,eta=.12)
        r=at_energy(H,vx,vy,.7,cfg)
        hp=2.5e-5
        gp=solve_self_energy(H,.7+hp,cfg,r['sigma'])['G']
        gm=solve_self_energy(H,.7-hp,cfg,r['sigma'])['G']
        self.assertLess(np.max(abs((gp-gm)/(2*hp)-r['derivative'])),8e-6)

    def test_clean_zero_W_matches_v38_primitive(self):
        d=sea_energy_scan(mesh=17,W=0,mu=0.,eta=.12,quadrature=240)
        self.assertLess(abs(d['dressed_sea']-d['clean_reference_if_W0']),1e-4)
        self.assertAlmostEqual(d['derivative_correction'],0.,delta=1e-12)

    def test_energy_integral_convergence_gapped(self):
        a=sea_energy_scan(mesh=17,W=.12,mu=0,quadrature=192)
        b=sea_energy_scan(mesh=17,W=.12,mu=0,quadrature=256)
        self.assertLess(abs(a['dressed_sea']-b['dressed_sea']),1e-4)
        self.assertLess(abs(b['dressed_sea']+0.841517),5e-4)
        self.assertLess(b['derivative_correction'],-.02)

    def test_energy_integral_convergence_metal(self):
        a=sea_energy_scan(mesh=17,W=.12,mu=1.1,quadrature=192)
        b=sea_energy_scan(mesh=17,W=.12,mu=1.1,quadrature=256)
        self.assertLess(abs(a['dressed_sea']-b['dressed_sea']),2e-4)
        self.assertGreater(abs(b['derivative_correction']),.1)
        self.assertEqual(b['protocol'],'energy_dependent_SCBA_bare_current_antisymmetric_sea_ONLY')

    def test_exact_disorder_W_zero_equals_single_Bloch_mesh(self):
        L=4
        k=-np.pi+2*np.pi*np.arange(L)/L
        x,y=np.meshgrid(k,k,indexing='ij')
        d=np.stack([np.sin(x),np.sin(y),-1+np.cos(x)+np.cos(y)])
        dx=np.stack([np.cos(x),np.zeros_like(x),-np.sin(x)])
        dy=np.stack([np.zeros_like(x),np.cos(y),-np.sin(y)])
        H=np.einsum('a...,aij->...ij',d,PAULI)
        vx=np.einsum('a...,aij->...ij',dx,PAULI)
        vy=np.einsum('a...,aij->...ij',dy,PAULI)
        target=2*np.pi*np.mean(sea_per_k(H,vx,vy,mu=0,eta=.12))
        exact=exact_finite_disorder_sea(L=4,W=0,mu=0,eta=.12,seeds=3)
        self.assertAlmostEqual(exact['mean'],target,delta=1e-11)
        self.assertAlmostEqual(exact['standard_deviation'],0.,delta=1e-12)

    def test_exact_disorder_reproducibility_and_finite_ensemble(self):
        a=exact_finite_disorder_sea(L=4,W=.12,seeds=6)
        b=exact_finite_disorder_sea(L=4,W=.12,seeds=6)
        self.assertEqual(a['samples'],b['samples'])
        self.assertGreater(a['standard_deviation'],.005)
        self.assertLess(a['standard_deviation'],.2)
        self.assertLess(abs(a['mean']-(-1.04)),.2)

    def test_disorder_sample_not_equal_scba_claim(self):
        ex=exact_finite_disorder_sea(L=4,W=.12,seeds=6)
        scba=sea_energy_scan(mesh=17,W=.12,mu=0,quadrature=224)
        self.assertNotEqual(ex['protocol'],scba['protocol'])
        self.assertIn('transverse current ladder',scba['explicit_missing'])
        self.assertGreater(abs(ex['mean']-scba['dressed_sea']),.05)

    def test_invalid_inputs_fail_closed(self):
        with self.assertRaises(ValueError):SCBAConfig(eta=0)
        with self.assertRaises(ValueError):SCBAConfig(mass=0)
        with self.assertRaises(ValueError):SCBAConfig(disorder_variance=-.1)
        with self.assertRaises(ValueError):sea_energy_scan(lower=1.,mu=0.)
        with self.assertRaises(ValueError):sea_energy_scan(quadrature=0)
        with self.assertRaises(ValueError):exact_finite_disorder_sea(L=2)
        with self.assertRaises(ValueError):exact_finite_disorder_sea(W=-1.)

if __name__=='__main__':unittest.main()
