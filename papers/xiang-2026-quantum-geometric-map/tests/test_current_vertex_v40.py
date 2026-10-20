"""Independent v40 model acceptance; this is NOT full disorder Hall conductivity."""
import unittest
import numpy as np
from qgmap.streda_clean_sea_v38 import qwz_model,I,SX
from qgmap.scba_streda_v39 import SCBAConfig,solve_self_energy
from qgmap.current_vertex_v40 import current_audit,pauli_vertex,exact_disorder_ra_torus

class V40CurrentVertex(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metal=current_audit(mesh=21,W=.12,energy=1.1)
        cls.gap=current_audit(mesh=21,W=.12,energy=0.)

    def test_scba_and_all_channel_residuals(self):
        self.assertLess(self.metal['selfenergy_residual'],1e-10)
        self.assertLess(self.metal['max_bethe_salpeter_residual'],1e-11)
        self.assertLess(self.metal['max_linear_condition'],500)
        self.assertEqual(set(self.metal['channel_bubbles']),{'ra','ar','rr','aa'})

    def test_density_ra_ward_matched_scba(self):
        self.assertLess(self.metal['ra_density_ward_residual'],1e-9)
        self.assertGreater(self.metal['ra_density_vertex_norm'],1)

    def test_ra_current_ladder_nonzero_metal(self):
        self.assertGreater(self.metal['current_correction_norms']['ra']['x'],.1)
        self.assertGreater(abs(self.metal['ra_vertex_effect']),.1)
        self.assertGreater(abs(self.metal['ra_antisymmetric_dressed']),1.)
        self.assertAlmostEqual(self.metal['ra_vertex_effect'],-.15220610119409872,delta=1e-7)

    def test_gap_current_ladder_nonzero(self):
        self.assertGreater(self.gap['current_correction_norms']['ra']['x'],.005)
        self.assertLess(abs(self.gap['ra_vertex_effect']),.001)

    def test_full_spectral_hall_antisym_zero(self):
        for a in (self.metal,self.gap):
            self.assertLess(a['four_channel_hall_antisym_abs'],1e-11)
            self.assertLess(a['bare_spectral_hall_antisym_abs'],1e-12)
            self.assertAlmostEqual(a['four_channel_spectral_xy']['real'],a['four_channel_spectral_yx']['real'],delta=1e-11)

    def test_cross_channel_hermitian_symmetric(self):
        c=self.metal['channel_bubbles']
        self.assertAlmostEqual(c['ra']['xy']['real'],c['ar']['yx']['real'],delta=1e-10)
        self.assertAlmostEqual(c['rr']['xy']['real'],c['rr']['yx']['real'],delta=1e-10)
        self.assertAlmostEqual(c['aa']['xy']['real'],c['aa']['yx']['real'],delta=1e-10)

    def test_rr_current_bz_ward_convergence(self):
        # A momentum-independent scalar SCBA Sigma implies <GR v GR>=<d_k GR>=0.
        # A finite shifted grid can have exponentially small sampling residuals.
        errors=[]
        for n in (21,31,41,61):
            a=current_audit(mesh=n,W=.12,energy=1.1)
            errors.append(max(a['current_correction_norms']['rr'].values()))
        self.assertGreater(errors[0],1e-4)
        self.assertLess(errors[-1],1e-8)
        self.assertTrue(all(a>b for a,b in zip(errors,errors[1:])))

    def test_clean_W_zero_all_current_ladders_vanish(self):
        a=current_audit(mesh=21,W=0.,energy=.7)
        for ch in a['current_correction_norms'].values():
            self.assertEqual(ch['x'],0.)
            self.assertEqual(ch['y'],0.)
        self.assertEqual(a['ra_vertex_effect'],0.)

    def test_su2_basis_covariance(self):
        H,vx,vy=qwz_model(17,-1)
        s=solve_self_energy(H,1.1,SCBAConfig(mesh=17,disorder_variance=.12))
        gr=s['G'];ga=gr.conj().swapaxes(-1,-2)
        a=pauli_vertex(gr,ga,vy,.12)
        theta=.39
        U=np.cos(theta)*I+1j*np.sin(theta)*SX
        uc=U.conj().T
        r=lambda x:uc@x@U
        b=pauli_vertex(r(gr),r(ga),r(vy),.12)
        np.testing.assert_allclose(b['constant'],r(a['constant']),atol=2e-11)
        self.assertLess(a['residual'],1e-11)
        self.assertLess(b['residual'],1e-11)

    def test_exact_torus_clean_matches_bloch_ra(self):
        # W=0 exact L=4 torus spectrum is the unshifted QWZ BZ mesh.
        e=exact_disorder_ra_torus(L=4,W=0.,seeds=3)
        self.assertAlmostEqual(e['mean'],-1.9355212426658206,delta=1e-11)
        self.assertEqual(e['standard_deviation'],0.)

    def test_exact_torus_disorder_reproducible_not_SCBA(self):
        a=exact_disorder_ra_torus(L=4,W=.12,seeds=6)
        b=exact_disorder_ra_torus(L=4,W=.12,seeds=6)
        self.assertEqual(a['samples'],b['samples'])
        self.assertGreater(a['standard_deviation'],.01)
        self.assertNotEqual(a['protocol'],self.metal['protocol'])
        self.assertIn('not SCBA',a['caveat'])

    def test_bse_fail_closed(self):
        with self.assertRaises(ValueError):pauli_vertex(np.ones((3,3,2,2)),np.ones((3,3,2,2)),np.ones((2,2)),.1)
        with self.assertRaises(ValueError):pauli_vertex(np.zeros((3,3,2,2)),np.zeros((3,3,2,2)),np.zeros((3,3,2,2)),-1.)
        with self.assertRaises(ValueError):exact_disorder_ra_torus(L=2)
        with self.assertRaises(ValueError):current_audit(mesh=8)

if __name__=='__main__':unittest.main()
