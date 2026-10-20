import unittest
import numpy as np
from qgmap import dirac, rashba, liouville, symmetry, orbital

class CoreTests(unittest.TestCase):
    def test_dirac_dispersion(self):
        p=dirac.Parameters(); kx=.03;ky=-.02
        ref=dirac.energies(kx,ky,p)
        for h in ("main-eq10","sm-eq62"):
            np.testing.assert_allclose(np.linalg.eigvalsh(dirac.hamiltonian(kx,ky,p,h)),ref,atol=1e-13)
    def test_zeeman_geometry_and_gauge(self):
        k=(.023,.017)
        for h in ("main-eq10","sm-eq62"):
            ref=dirac.local_zbc_analytic(*k,h)
            got=dirac.local_zbc_eigensolver(*k,helicity=h,phases=(1.0,2.7))
            np.testing.assert_allclose(got,ref,atol=1e-12,rtol=1e-12)
    def test_opposite_helicity(self):
        a=np.array(dirac.local_zbc_eigensolver(.023,.017,helicity="main-eq10"))
        b=np.array(dirac.local_zbc_eigensolver(.023,.017,helicity="sm-eq62"))
        np.testing.assert_allclose(a,-b,atol=1e-12)
    def test_formula_signed_branches(self):
        for beta in (.03,.125,.35):
            a=float(dirac.eq11_kernel(beta));b,err=dirac.eq63_literal_integral(beta)
            self.assertLess(abs(a+b),1e-10);self.assertLess(err,1e-11)
    def test_source_beta(self):
        self.assertAlmostEqual(float(dirac.eq11_kernel(1/8)),-.0025064159066917205,places=13)
    def test_tanh_thermal_step(self):
        m=np.array([-.01,0,.01]);t=20
        expected=dirac.eq11_kernel(.125)*np.tanh(m/(2*dirac.KB_EV_K*t))
        np.testing.assert_allclose(dirac.continuum_temperature_kernel(m,t,.125),expected)
    def test_shift_control(self):
        self.assertLess(abs(dirac.equilibrium_drude_shift_invariant()),1e-12)
    def test_zeeman_prefactor(self):
        self.assertAlmostEqual(dirac.zeeman_time_prefactor(),.4397002,places=5)
    def test_rashba_two_rings(self):
        m=1.;lam=.4;mu=-.04
        a=rashba.zero_temperature(m,lam,mu)
        b=rashba.direct_fermi_ring_sum(m,lam,mu)
        for key in ("bare","metric","total"):
            self.assertLess(abs(a[key]-b[key]),1e-12)
        self.assertAlmostEqual(a["total"],0,places=14)
    def test_rashba_positive(self):
        a=rashba.zero_temperature(1,.4,.5)
        b=rashba.direct_fermi_ring_sum(1,.4,.5)
        for key in ("bare","metric","total"):
            self.assertLess(abs(a[key]-b[key]),1e-12)
    def test_parabolic_ohe(self):
        for m in (.5,1.,2.):
            out=orbital.ohe_parabolic_2d(m,.2,300)
            self.assertLess(out.relative_error,2e-10)
            self.assertAlmostEqual(out.sigma_xy_over_tau2,-out.sigma_yx_over_tau2,places=12)
    def test_symmetry_table(self):
        rows=symmetry.table_i_acceptance()
        self.assertEqual(len(rows),11)
        self.assertTrue(all(v["passed"] for v in rows.values()))
    def test_liouville_source_resolvent(self):
        h=np.diag([-1.,.4]).astype(complex)
        v=np.array([[.1,.17],[.17,-.2]],complex)
        b=np.array([[.2,.13j],[-.13j,-.1]],complex)
        rho=np.diag([.9,.2]).astype(complex)
        out=liouville.two_color_response(h,v,b,rho,omega_e=.24,omega_b=.13,eta_e=.04,eta_b=.03)
        def liouv(x):return h@x-x@h
        source_e=v@rho-rho@v;source_b=b@rho-rho@b
        self.assertLess(np.max(abs((.24+.04j)*out["rho_E"]-liouv(out["rho_E"])-source_e)),1e-13)
        self.assertLess(np.max(abs((.13+.03j)*out["rho_B"]-liouv(out["rho_B"])-source_b)),1e-13)
        source=v@out["rho_B"]-out["rho_B"]@v+b@out["rho_E"]-out["rho_E"]@b
        self.assertLess(np.max(abs((.37+.07j)*out["total"]-liouv(out["total"])-source)),1e-12)
        np.testing.assert_allclose(out["total"],out["B_then_E"]+out["E_then_B"])
    def test_population_separate(self):
        x=liouville.equilibrium_population_derivative(np.array([-.1,.2]),np.array([.3,-.4]),0,.05)
        self.assertEqual(x.shape,(2,))
        self.assertLess(x[0],0);self.assertGreater(x[1],0)
    def test_unsafe_parameter_rejected(self):
        with self.assertRaises(ValueError):dirac.eq11_kernel(1.)
        with self.assertRaises(ValueError):dirac.local_zbc_eigensolver(0,0)
        with self.assertRaises(ValueError):rashba.zero_temperature(1,.4,0)
        with self.assertRaises(ValueError):liouville.two_color_response(np.eye(2),np.eye(2),np.eye(2),np.eye(2),eta_e=0)
if __name__=="__main__":unittest.main()
