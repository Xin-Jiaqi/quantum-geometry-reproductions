import unittest
import numpy as np
from igmc_lab.core import dirac_hamiltonian,dirac_grid,pair_kernels,occupation
from igmc_lab.green_wavepacket_audit import green_eq44_resolvent,wavepacket_cos_sin,analytical_cosine_harmonics

class Round8GreenWavepacket(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(839)
        self.kx=rng.uniform(-2.4,2.4,14)
        self.ky=rng.uniform(-2.4,2.4,14)
        self.H,self.dH=dirac_hamiltonian(self.kx,self.ky)
    def test_lesser_green_vs_original_kubo_independent(self):
        for omega in (.003,.019,.071):
            for eta in (0.,.00013):
                for axis,spin in [('x','x'),('x','y'),('y','z')]:
                    observed=green_eq44_resolvent(self.H,self.dH[axis],spin,omega,eta=eta)
                    eig,U=np.linalg.eigh(self.H)
                    V=U.conj().swapaxes(-2,-1)@self.dH[axis]@U
                    from igmc_lab.core import SPINS
                    S=U.conj().swapaxes(-2,-1)@SPINS[spin]@U
                    a=V[...,0,1]*S[...,1,0]
                    delta=eig[...,1]-eig[...,0]
                    expected=a/(omega-delta+1j*eta)-a.conj()/(omega+delta+1j*eta)
                    self.assertLess(float(np.max(abs(observed-expected))),4e-13)
    def test_wavepacket_vs_lesser_green_all_quadratures(self):
        for omega in (.003,.025,.07):
            for ax,sp in [('x','x'),('x','y'),('y','x'),('y','z')]:
                wave=wavepacket_cos_sin(self.H,self.dH[ax],sp,omega)
                green=green_eq44_resolvent(self.H,self.dH[ax],sp,omega)
                self.assertLess(float(np.max(abs(wave['K_from_wavepacket']-green))),1e-12)
                self.assertLess(float(np.max(abs(wave['K_main_lowfrequency']-green))),.2)
    def test_sm_eq34_metric_q_is_double_of_exact_wavepacket(self):
        result=wavepacket_cos_sin(self.H,self.dH['x'],'x',.007)
        gap=np.diff(np.linalg.eigvalsh(self.H),axis=-1)[...,0]
        exact=result['P_cos_per_b']
        # Correct strict omega->0 amplitude is +2Q/Delta; printed SM intermediate
        # uses Qsm=2ReA, so yields twice the correct cosine amplitude.
        leading=2*result['Q_main']/gap
        self.assertLess(float(np.max(abs(exact-leading))),2e-3)
        self.assertLess(float(np.max(abs(2*leading-2*result['Q_main']*2/gap))),1e-14)
        self.assertGreater(float(np.max(abs(2*leading-exact))),.05)
    def test_main_eq10_sine_sign_vs_exact_wavepacket(self):
        wave=wavepacket_cos_sin(self.H,self.dH['x'],'x',.02)
        assert np.max(abs(wave['Z_main']))>1e-3
        # Main Eq(10) F=-omega*Z/Delta^2 with P=-b*F*sin => Psin=+omega*Z/Delta^2.
        # SM eq(24) wavepacket instead gives Psin=-omega*Z/(Delta^2-omega^2).
        gap=np.diff(np.linalg.eigvalsh(self.H),axis=-1)[...,0]
        calculated=wave['P_sin_per_b']
        main=wave['P_sin_main_Eq10_per_b']
        self.assertLess(float(np.max(abs(calculated+main/(1-(.02/gap)**2)))),2e-13)
        self.assertGreater(float(np.max(abs(calculated-main))),1e-3)
    def test_cosine_fourier_coefficient_is_half(self):
        a,b=analytical_cosine_harmonics()
        self.assertAlmostEqual(a,.5,places=14)
        self.assertAlmostEqual(b,.5,places=14)
    def test_wavepacket_not_valid_for_metal(self):
        with self.assertRaises(ValueError):wavepacket_cos_sin(self.H,self.dH['x'],'x',.003,mu=-10.)
    def test_green_no_intrinsic_zero_frequency(self):
        with self.assertRaises(ValueError):green_eq44_resolvent(self.H,self.dH['x'],'y',0.)
if __name__=='__main__':unittest.main()
