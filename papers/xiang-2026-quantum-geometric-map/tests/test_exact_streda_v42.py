"""v42 microscopic disordered Streda gates independent of v41 Lehmann definitions."""
import unittest
import numpy as np
from qgmap.exact_hall_v41 import torus,hall_lehmann
from qgmap.exact_streda_v42 import (spectral_streda,direct_surface_trace,
    direct_sea_quadrature,streda_audit,ensemble_audit)

class ExactFiniteDisorderStreda(unittest.TestCase):
    def test_clean_and_disordered_gap_sum_to_full_lehmann(self):
        for W in (0.,.12):
            H,x,y=torus(4,-1.,W,1947,(.11,.33))
            target=hall_lehmann(H,x,y,mu=0.)['hall_real']
            errors=[]
            for eta in (.16,.05,.01,.003):
                v=spectral_streda(H,x,y,0.,eta)
                errors.append(abs(v['total']-target))
            self.assertLess(errors[-1],2e-7)
            self.assertLess(errors[-1],errors[0]*1e-4)
    def test_metal_has_large_opposite_surface_and_sea(self):
        H,x,y=torus(4,-1.,.12,1947,(.11,.33))
        finite=spectral_streda(H,x,y,mu=1.1,eta=.05)
        self.assertLess(finite['surface_Ia'],-1.)
        self.assertGreater(finite['sea_II'],.5)
        self.assertLess(abs(finite['total']),.5)
        target=hall_lehmann(H,x,y,mu=1.1)['hall_real']
        small=spectral_streda(H,x,y,mu=1.1,eta=.001)
        self.assertLess(abs(small['total']-target),2e-6)
        self.assertGreater(abs(finite['sea_II']-target),.5)
    def test_matrix_trace_surface_independent_of_pair_formula(self):
        for W in (0.,.12):
            H,x,y=torus(3,-1.,W,1947,(.11,.33))
            sp=spectral_streda(H,x,y,1.1,.09)
            direct=direct_surface_trace(H,x,y,1.1,.09)
            self.assertAlmostEqual(sp['surface_Ia'],direct['surface_Ia_hall_odd'],delta=2e-11)
            self.assertLess(direct['surface_Ib_hall_odd_complex_abs'],1e-10)
    def test_matrix_energy_quad_independent_of_primitive(self):
        for mu in (0.,1.1):
            H,x,y=torus(3,-1.,.12,1947,(.11,.33))
            analytical=spectral_streda(H,x,y,mu,.14)['sea_II']
            for lower in (-30.,-55.):
                quad=direct_sea_quadrature(H,x,y,mu,.14,lower)
                self.assertLess(abs(quad['sea_II']-analytical),3e-7)
    def test_v38_clean_bloch_sea_matches_exact_realspace_v42(self):
        from qgmap.streda_clean_sea_v38 import qwz_model,sea_per_k
        L=5;theta=(.11,.33);offset=tuple(t/(2*np.pi) for t in theta)
        H,x,y=torus(L,-1.,0.,1947,theta)
        h,vx,vy=qwz_model(L,-1.,offset=offset)
        for mu in (0.,1.1):
            for eta in (.09,.025):
                exact=spectral_streda(H,x,y,mu,eta)['sea_II']
                clean=float(2*np.pi*np.mean(sea_per_k(h,vx,vy,mu,eta)))
                self.assertAlmostEqual(exact,clean,delta=5e-11)

    def test_xy_antisymmetry_and_constant_unitary_basis(self):
        H,x,y=torus(3,-1.,.12,1997,(.17,.23))
        original=spectral_streda(H,x,y,1.1,.07)
        flip=spectral_streda(H,y,x,1.1,.07)
        self.assertAlmostEqual(flip['total'],-original['total'],delta=2e-11)
        rng=np.random.default_rng(41)
        z=rng.normal(size=H.shape)+1j*rng.normal(size=H.shape)
        q,r=np.linalg.qr(z)
        phase=np.diag(r);u=q*np.where(abs(phase)>0,phase/np.maximum(abs(phase),1e-30),1)[None,:]
        transform=lambda a:u.conj().T@a@u
        got=spectral_streda(transform(H),transform(x),transform(y),1.1,.07)
        self.assertAlmostEqual(got['surface_Ia'],original['surface_Ia'],delta=1e-10)
        self.assertAlmostEqual(got['sea_II'],original['sea_II'],delta=1e-10)
    def test_exterior_phase_chern_signs(self):
        for m,sgn in ((-1,-1),(1,1),(-3,0),(3,0)):
            H,x,y=torus(3,m,0.,17,(.11,.33))
            a=spectral_streda(H,x,y,0.,.002)
            l=hall_lehmann(H,x,y,mu=0.)['hall_real']
            self.assertLess(abs(a['total']-l),1e-5)
            if sgn!=0:self.assertGreater(sgn*l,0)
    def test_exact_finite_seed_average_not_scba(self):
        r=ensemble_audit(L=3,W=.12,mu=1.1,eta=.002,
                         seeds=(1947,1948,1949,1950))
        self.assertEqual(r['nsamples'],4)
        self.assertAlmostEqual(r['mean_surface_Ia']+r['mean_sea_II'],
                               r['mean_streda_total'],delta=1e-12)
        self.assertLess(abs(r['mean_total_minus_lehmann']),.002)
        self.assertEqual(r['scientific_status'],'PARTIAL')
    def test_audit_protocol_and_bounds(self):
        v=streda_audit(L=3,W=.12,mu=0.,etas=(.10,.025,.006))
        self.assertEqual(v['scientific_status'],'PARTIAL')
        self.assertLess(v['rows'][-1]['abs_error_vs_exact_lehmann'],3e-6)
        self.assertGreater(v['nearest_level_distance'],.1)
    def test_fail_closed(self):
        H,x,y=torus(3,-1.,.12,1947,(.11,.33))
        for bad in (0.,-.01,np.nan):
            with self.assertRaises(ValueError):spectral_streda(H,x,y,0.,bad)
        with self.assertRaises(ValueError):spectral_streda(H,x,np.eye(2),0.,.05)
        with self.assertRaises(ValueError):direct_sea_quadrature(H,x,y,0.,.05,lower=1.)
        with self.assertRaises(ValueError):ensemble_audit(seeds=(1,1))

if __name__=='__main__':unittest.main()
