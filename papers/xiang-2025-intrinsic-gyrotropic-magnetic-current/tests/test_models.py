import unittest
import numpy as np
from igmc_lab.core import (SX,SY,SZ,dirac_grid,dirac_hamiltonian,
      dirac_kernel,dirac_response,fermi_surface_response,rashba_kernel,
      rashba_analytic,rashba_integrate,kernel_response,pair_kernels)

class PaperModelTests(unittest.TestCase):
    def test_dirac_gamma_gap(self):
        H,_=dirac_hamiltonian(np.array([0.]),np.array([0.]))
        self.assertAlmostEqual(np.linalg.eigvalsh(H)[0,0],-.25)
        self.assertAlmostEqual(np.linalg.eigvalsh(H)[0,1],+.25)

    def test_rashba_band_minimum_and_flat_top(self):
        alpha=3.80998212/1.1;lam=.6
        self.assertAlmostEqual(-lam**2/(4*alpha),-.025985,delta=5e-6)
        self.assertAlmostEqual(rashba_analytic(.001),lam/(4*np.pi*alpha))
        self.assertEqual(rashba_analytic(-.05),0.)

    def test_rashba_analytic_vs_GEOMETRIC_integral(self):
        for mu in [-.015,0.,.009]:
            observed=rashba_integrate(mu,nrad=1200,nphi=64)
            expect=rashba_analytic(mu)
            self.assertAlmostEqual(observed,expect,delta=.0002)

    def test_rashba_analytic_angular_texture(self):
        k=.17
        # For n=lower, m=upper: Z_xy ~ k_y^2/k^3, Z_xx ~ +/-kx*ky/k^3.
        p=rashba_kernel(np.array([0.,k,k/2]),np.array([k,0.,k/2]))
        self.assertAlmostEqual(p['Z01'][0],1/k,places=10)
        self.assertAlmostEqual(p['Z01'][1],0,places=10)
        self.assertAlmostEqual(p['Z01'][2],0.5/np.sqrt(2*(k/2)**2),places=10)
        self.assertAlmostEqual(p['Q01'][0],0,places=10)

    def test_gauge_covariance_eigenvectors_at_fixed_k(self):
        kx=np.array([.13]);ky=np.array([.20])
        H,dH=dirac_hamiltonian(kx,ky)
        k0=pair_kernels(H,dH,'x','y')
        rng=np.random.default_rng(2);A=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        U,_=np.linalg.qr(A)
        rot=lambda x: U.conj().T@x@U
        H2=rot(H);dH2={a:rot(h) for a,h in dH.items()}
        # Spin matrices transform too; standard pair_kernels uses fixed physical Pauli,
        # so test the invariant product directly with unitary conjugated operators.
        e,V=np.linalg.eigh(H[0]);ev0=V.conj().T@dH['x'][0]@V
        evs0=V.conj().T@SY@V
        e,W=np.linalg.eigh(H2[0]);ev1=W.conj().T@dH2['x'][0]@W
        evs1=W.conj().T@rot(SY)@W
        self.assertAlmostEqual(abs(1j*ev0[0,1]*evs0[1,0]/(e[1]-e[0])-1j*ev1[0,1]*evs1[1,0]/(e[1]-e[0])),0,places=12)
        self.assertTrue(np.isfinite(k0['Z01']).all())

    def test_dirac_insulator_conduction_zero_displacement_finite(self):
        C=dirac_response(51,mu=0.,alpha='x',beta='y')
        D=dirac_response(51,mu=0.,alpha='x',beta='x')
        self.assertLess(abs(C['sigma_C']),1e-7)  # convergence tracked separately
        C101=dirac_response(101,mu=0.,alpha='x',beta='y')
        self.assertLess(abs(C101['sigma_C']),1e-11)
        self.assertGreater(abs(D['sigma_D_per_omega']),1e-4)
        self.assertLess(abs(D['Kubo_real']-D['sigma_C']),1e-8)
        self.assertAlmostEqual(D['Kubo_imag'],D['sigma_D'],delta=2e-7)

    def test_dirac_C_fermi_surface_equivalence(self):
        kx,ky=dirac_grid(121)
        pair=dirac_kernel(kx,ky,'x','y')
        for mu in [.53,.85]:
            sf=fermi_surface_response(pair,mu,temperature_K=200)
            geom=kernel_response(pair,mu,temperature_K=200,omega_eV=.001)['sigma_C']
            self.assertAlmostEqual(sf,geom,delta=2e-5)

    def test_antiunitary_rashba_C_nonzero_D_cancel(self):
        theta=np.linspace(.1,2*np.pi+.1,64,endpoint=False)
        kx=.2*np.cos(theta);ky=.2*np.sin(theta)
        C=rashba_kernel(kx,ky,'x','y')
        D=rashba_kernel(kx,ky,'x','x')
        self.assertGreater(abs(np.mean(C['Z01'])),.2)
        self.assertLess(abs(np.mean(D['Q01'])),1e-12)

    def test_rashba_explicit_time_reversal(self):
        from igmc_lab.core import rashba_hamiltonian
        kx=np.array([.113,-.097]);ky=np.array([-.054,.082])
        H,_=rashba_hamiltonian(kx,ky)
        Hminus,_=rashba_hamiltonian(-kx,-ky)
        transformed=SY@H.conj()@SY
        self.assertLess(float(np.max(abs(Hminus-transformed))),1e-14)

    def test_frequency_slope(self):
        r1=dirac_response(51,mu=0.,alpha='x',beta='x',omega_eV=.0008)
        r2=dirac_response(51,mu=0.,alpha='x',beta='x',omega_eV=.0016)
        self.assertAlmostEqual(r2['sigma_D'],2*r1['sigma_D'],places=14)
        self.assertAlmostEqual(r1['Kubo_imag']/r1['sigma_D'],1,delta=2e-5)
        self.assertAlmostEqual(r2['Kubo_imag']/r2['sigma_D'],1,delta=4e-5)

    def test_gap_fail_closed(self):
        kx=np.array([0.]);ky=np.array([0.])
        H,dH=dirac_hamiltonian(kx,ky,beta=0.,M=0.)
        with self.assertRaises(ValueError):pair_kernels(H,dH,'x','y')
        with self.assertRaises(ValueError):dirac_response(31,omega_eV=2.)

if __name__=='__main__':unittest.main()
