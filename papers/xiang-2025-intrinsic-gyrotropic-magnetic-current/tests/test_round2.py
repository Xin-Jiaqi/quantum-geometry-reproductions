"""Round 2 physical-contract tests; separate from frozen round-1 data."""
import math
import unittest
import numpy as np
from igmc_lab.core import (
    SY, dirac_grid, dirac_hamiltonian, dirac_kernel, dirac_response,
    rashba_analytic, energy_from_frequency, physical_response_si,
)

class TimeReversalAndPhysicalUnits(unittest.TestCase):
    def test_dirac_time_reversed_hamiltonian_exact(self):
        kx=np.array([.21,-.09,.53]);ky=np.array([-.32,.12,.35])
        H,_=dirac_hamiltonian(kx,ky,mass_sign=1)
        HT,_=dirac_hamiltonian(-kx,-ky,mass_sign=-1)
        self.assertLess(float(np.max(abs(HT-SY@H.conj()@SY))),1e-14)

    def test_local_z_even_q_odd_under_tr(self):
        kx,ky=dirac_grid(31)
        for alpha,beta in [('x','y'),('x','x')]:
            p=dirac_kernel(kx,ky,alpha,beta,mass_sign=1)
            t=dirac_kernel(-kx,-ky,alpha,beta,mass_sign=-1)
            self.assertLess(np.max(abs(p['Z01']-t['Z01'])),1e-12)
            self.assertLess(np.max(abs(p['Q01']+t['Q01'])),1e-12)

    def test_conduction_even_displacement_odd_under_tr(self):
        for mu in [0.,.85]:
            xp=dirac_response(71,mu,alpha='x',beta='y',mass_sign=1)
            xm=dirac_response(71,mu,alpha='x',beta='y',mass_sign=-1)
            dp=dirac_response(71,mu,alpha='x',beta='x',mass_sign=1)
            dm=dirac_response(71,mu,alpha='x',beta='x',mass_sign=-1)
            self.assertAlmostEqual(xp['sigma_C'],xm['sigma_C'],delta=1e-12)
            self.assertAlmostEqual(dp['sigma_D_per_omega'],-dm['sigma_D_per_omega'],delta=1e-12)
        # In the metal this C is not trivially zero; in the gap D not trivially zero.
        self.assertGreater(abs(xp['sigma_C']),.01)
        self.assertGreater(abs(dirac_response(71,0,alpha='x',beta='x')['sigma_D_per_omega']),.01)

    def test_kubo_complex_partner_parity(self):
        rp=dirac_response(81,.85,alpha='x',beta='y',mass_sign=1)
        rm=dirac_response(81,.85,alpha='x',beta='y',mass_sign=-1)
        dp=dirac_response(81,0.,alpha='x',beta='x',mass_sign=1)
        dm=dirac_response(81,0.,alpha='x',beta='x',mass_sign=-1)
        self.assertAlmostEqual(rp['Kubo_real'],rm['Kubo_real'],delta=1e-12)
        self.assertAlmostEqual(dp['Kubo_imag'],-dm['Kubo_imag'],delta=1e-12)

    def test_thz_frequency_defined_by_sm_section_vi(self):
        E=energy_from_frequency(1e12,protocol='cyclic_hz')
        self.assertAlmostEqual(E,.0041356676944,delta=2e-11)
        self.assertAlmostEqual(E/energy_from_frequency(1e12,protocol='angular_rad_s'),2*math.pi,places=12)
        with self.assertRaises(ValueError):energy_from_frequency(1e12,protocol='Hz')
        with self.assertRaises(ValueError):energy_from_frequency(-1,protocol='cyclic_hz')

    def test_si_prefactor_and_rashba_absolute_scale(self):
        # Eq.16 versus main Figure 1(b), without fit parameters, for g=2.
        self.assertAlmostEqual(physical_response_si(1.),281.79403263785,delta=1e-7)
        self.assertAlmostEqual(physical_response_si(rashba_analytic(0.)),3.88456976939,delta=1e-8)
        self.assertAlmostEqual(physical_response_si(.013,g_factor=4),2*physical_response_si(.013,g_factor=2),places=10)
        with self.assertRaises(ValueError):physical_response_si(1,g_factor=-1)

    def test_direct_eq16_in_si_independent_of_k_space_integral(self):
        from igmc_lab.core import SI_E_C,SI_HBAR_J_S,SI_MU_B_J_T
        lam_J_meter=.6*SI_E_C*1e-10
        effective_mass_kg=1.1*9.1093837139e-31
        direct_formula = (SI_E_C*2*SI_MU_B_J_T*lam_J_meter*effective_mass_kg
                          /(2*math.pi*SI_HBAR_J_S**3))
        self.assertAlmostEqual(direct_formula,physical_response_si(rashba_analytic(0.)),delta=2e-8)

    def test_dirac_displacement_absolutely_not_paper_matched(self):
        energy=energy_from_frequency(1e12,protocol='cyclic_hz')
        d=dirac_response(101,0.,omega_eV=energy,alpha='x',beta='x')
        self.assertAlmostEqual(physical_response_si(d['sigma_D']),.01775389322,delta=2e-9)
        # No arbitrary scaling to the article's Fig. 2(d) ~0.04 A/(m T) plateau.

if __name__=='__main__':unittest.main()
