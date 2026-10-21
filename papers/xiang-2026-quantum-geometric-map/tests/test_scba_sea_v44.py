"""Bounded v44 SCBA sea checks: not full conserving Hall nor second-order EB."""
import unittest
import numpy as np
from qgmap.scba_sea_v44 import (sea_scan,sea_kernel,_setup,candidate_total,
    clean_bloch_pair,exact_ensemble_total,compare_candidate)
from qgmap.scba_streda_v39 import solve_self_energy
from qgmap.exact_streda_v42 import spectral_streda
from qgmap.exact_hall_v41 import torus


class SeaV44(unittest.TestCase):
    def test_clean_two_band_independent_pole_vs_integration(self):
        for n in (21,35):
            clean=clean_bloch_pair(n,mu=1.1,eta=.12)
            a=candidate_total(mesh=n,W=0,mu=1.1,eta=.12,
                              lower=-12,quadrature=420)
            self.assertAlmostEqual(a['surface_Ia'],clean['surface_Ia'],delta=2e-12)
            self.assertAlmostEqual(a['sea_II_candidate'],clean['sea_II'],delta=2e-6)
            self.assertEqual(a['scientific_status'],'PARTIAL')

    def test_ward_energy_fd_three_energies_and_two_disorders(self):
        for W in (.03,.12):
            H,x,y,c=_setup(mesh=17,W=W)
            for E in (-2.1,-.75,.82):
                a=sea_kernel(E,H,x,y,c,finite_difference_h=2e-4)
                self.assertLess(a['dG_finite_difference_error'],2e-4)
                self.assertLess(a['sigma_finite_difference_error'],4e-7)
                self.assertLess(a['gamma0_residual'],1e-9)
                self.assertLess(a['scba_residual'],1e-9)

    def test_self_energy_derivative_is_not_naively_identity(self):
        H,x,y,c=_setup(mesh=21,W=.12)
        p=sea_kernel(.82,H,x,y,c)
        self.assertGreater(p['gamma0_minus_I'],.1)

    def test_direct_finite_difference_integrand_independent_of_BSE(self):
        H,vx,vy,c=_setup(mesh=17,W=.12)
        E=.82;h=1e-5
        a=sea_kernel(E,H,vx,vy,c)
        z=solve_self_energy(H,E,c)
        gp=solve_self_energy(H,E+h,c,initial=z['sigma'])['G']
        gm=solve_self_energy(H,E-h,c,initial=z['sigma'])['G']
        d=(gp-gm)/(2*h)
        g=z['G']; ga=g.conj().swapaxes(-1,-2); da=d.conj().swapaxes(-1,-2)
        tr=lambda v:np.trace(v,axis1=-2,axis2=-1)
        val=.5*np.real(np.mean(tr(vx@g@vy@d-vx@d@vy@g-vx@ga@vy@da+vx@da@vy@ga)))
        self.assertAlmostEqual(val,a['sea_density'],delta=1e-5)

    def test_finite_disorder_quadrature_and_tail_convergence(self):
        for W in (.03,.12):
            a=sea_scan(mesh=21,W=W,quadrature=360,lower=-9)
            b=sea_scan(mesh=21,W=W,quadrature=420,lower=-12)
            self.assertLess(abs(a['sea_II_candidate']-b['sea_II_candidate']),3e-5)

    def test_scba_candidate_does_not_claim_conserving_hall(self):
        r=candidate_total(mesh=21,W=.12,quadrature=260)
        self.assertIn('candidate_total',r)
        self.assertNotIn('hall_physical_dc',r)
        self.assertEqual(r['scientific_status'],'PARTIAL')
        self.assertTrue(any('conserving' in x for x in r['not_proven']))

    def test_exact_disordered_streda_per_seed_is_maintained(self):
        L=3; W=.12; mu=1.1; eta=.12; seed=1947
        h,x,y=torus(L,-1.,W,seed,(0.,0.))
        v=spectral_streda(h,x,y,mu,eta)
        r=exact_ensemble_total(L=L,W=W,mu=mu,eta=eta,seeds=(seed,),ntheta=2)
        self.assertAlmostEqual(r['full_total'],r['surface_Ia']+r['sea_II'],delta=1e-12)
        self.assertGreater(abs(v['total']-r['full_total']),1e-5)

    def test_weak_moderate_disorder_candidate_vs_exact_ensemble(self):
        for W in (.03,.12):
            out=compare_candidate(mesh=21,L=3,W=W,mu=1.1,eta=.12,
                                  seeds=tuple(range(1947,1955)),ntheta=3,
                                  quadrature=300,lower=-12)
            self.assertEqual(out['scientific_status'],'PARTIAL')
            self.assertLess(abs(out['total_discrepancy']),.4)
            self.assertGreater(out['finite_exact_reference']['seed_sem_total'],0)

    def test_clean_reference_is_explicitly_not_finite_disorder(self):
        v=compare_candidate(mesh=21,L=3,W=0,quadrature=420,lower=-12)
        self.assertTrue(v['finite_exact_reference']['not_disorder_average'])
        self.assertLess(abs(v['total_discrepancy']),2e-6)
        with self.assertRaises(ValueError):exact_ensemble_total(W=0)

    def test_invalid_inputs_rejected(self):
        for kwargs in ({'W':-.01},{'eta':0},{'mesh':8},{'mu':0.,'lower':1.},
                       {'quadrature':5}):
            with self.assertRaises(ValueError):sea_scan(**kwargs)
        H,x,y,c=_setup(17,W=.12)
        with self.assertRaises(ValueError):sea_kernel(.2,H,x,y,c,finite_difference_h=-1e-3)
        with self.assertRaises(ValueError):exact_ensemble_total(seeds=(1,1))

if __name__=='__main__':unittest.main()
