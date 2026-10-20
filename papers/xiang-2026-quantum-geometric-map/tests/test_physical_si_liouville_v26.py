"""v26 direct SI 2x2 Liouville source, angular measure and h/hbar audit."""
import unittest, math
from dataclasses import replace
from scipy import constants as c
from qgmap.si_fig2_v25 import FigureInputs,evaluate,geometry,independent_angular_geometry
from qgmap.physical_si_liouville_v26 import (SIConfig,eq20_physical,eq11_SI_magnitude,
                                            two_color_physical,frequency_hypotheses,kgrid)

class PhysicalSI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.p=SIConfig(nr=92,nphi=128)
    def test_01_eq20_matches_restored_eq11(self):
        p=self.p
        actual=eq20_physical(p).real/p.E_Vm
        expected=eq11_SI_magnitude(p)
        self.assertLess(abs(abs(actual)/expected-1),1e-5)
        self.assertLess(actual,0)
    def test_02_no_hidden_factor_2pi(self):
        a=abs(eq20_physical(self.p)/self.p.E_Vm)
        target=eq11_SI_magnitude(self.p)
        self.assertAlmostEqual(a/target,1.,delta=1e-5)
        self.assertGreater(abs(a/target-1/(2*math.pi)),.8)
    def test_03_hbar_frequency_matches_tau_phase(self):
        x=frequency_hypotheses(self.p)
        self.assertAlmostEqual(x['omega_Z_rad_s'],2*math.pi*x['f_Z_Hz'],delta=.01)
        self.assertAlmostEqual(x['tau_omega_Z_rad'],.4397050020926946,places=11)
    def test_04_angular_integral_has_correct_2pi_square_measure(self):
        self.assertAlmostEqual(geometry(self.p.beta),independent_angular_geometry(self.p.beta),places=12)
        self.assertNotAlmostEqual(geometry(self.p.beta),2*math.pi*independent_angular_geometry(self.p.beta),places=6)
    def test_05_vf_scale_invariance(self):
        p=self.p
        for vf in (.22,.48):
            q=eq20_physical(replace(p,vF_eVnm=vf)).real
            self.assertAlmostEqual(q/eq20_physical(p).real,1.,delta=5e-5)
    def test_06_electric_field_linearity(self):
        p=self.p
        j=eq20_physical(p)
        twice=eq20_physical(replace(p,E_Vm=2*p.E_Vm))
        self.assertAlmostEqual(abs(twice/j),2.,delta=2e-9)
    def test_07_zeeman_field_linearity(self):
        p=self.p
        j=eq20_physical(p)
        twice=eq20_physical(replace(p,B_T=2*p.B_T))
        self.assertAlmostEqual(abs(twice/j),2.,delta=2e-9)
    def test_08_g_factor_linearity(self):
        p=self.p
        j=eq20_physical(p)
        twice=eq20_physical(replace(p,g=2*p.g))
        self.assertAlmostEqual(abs(twice/j),2.,delta=2e-9)
    def test_09_tau_linearity_eq20(self):
        p=self.p
        j=eq20_physical(p)
        twice=eq20_physical(replace(p,tau_s=2*p.tau_s))
        self.assertAlmostEqual(abs(twice/j),2.,delta=2e-9)
    def test_10_sign_mu_inversion(self):
        p=self.p
        plus=eq20_physical(p).real
        minus=eq20_physical(replace(p,mu_eV=-p.mu_eV)).real
        self.assertLess(abs(plus+minus)/abs(plus),1e-5)
    def test_11_sign_beta_inversion(self):
        p=self.p
        plus=eq20_physical(p).real
        minus=eq20_physical(replace(p,beta=-p.beta)).real
        self.assertLess(abs(plus+minus)/abs(plus),1e-5)
    def test_12_resolvent_pole_is_tau_effective(self):
        p=replace(self.p,nr=92,nphi=128)
        # All terms, not merely rho(EB;5); nonleading terms remain at finite eta.
        for scale in (2.,1.,.5):
            j=two_color_physical(replace(p,eta_scale=scale))['sum_A_per_m'].real
            expected=eq20_physical(p).real/scale
            self.assertLess(abs(j/expected-1),.013)
    def test_13_eta_branch_no_implicit_h(self):
        p=self.p
        x=two_color_physical(p)
        self.assertAlmostEqual(x['eta_energy_J'],c.hbar/p.tau_s,delta=1e-28)
        self.assertAlmostEqual(x['natural_dimensionless_tau_Zeeman_phase'],frequency_hypotheses(p)['tau_omega_Z_rad'],places=12)
    def test_14_wavevector_measure(self):
        p=self.p
        _,w,vf=kgrid(p)
        cutoff=p.kmax_eV*c.e/vf
        self.assertAlmostEqual(w.sum(),cutoff**2/(4*math.pi),delta=1e-8*cutoff**2)
    def test_15_independent_fig2_voltage_consistency(self):
        s=evaluate(FigureInputs())
        # There is no new 2pi in restored numerical Dirac response.
        self.assertAlmostEqual(eq11_SI_magnitude(self.p),abs(evaluate(FigureInputs(),mu_mev=60.,temp_K=self.p.thermal_eV*c.e/c.Boltzmann)['sigma_sheet_S']),delta=1e-16)
        self.assertAlmostEqual(abs(s['open_circuit_voltage_microV']),268.262,delta=.002)
        self.assertAlmostEqual(abs(evaluate(FigureInputs(),hypothesis='single_2pi_diagnostic')['open_circuit_voltage_microV']),42.695,delta=.002)
    def test_16_mu_zero_zero_response(self):
        p=replace(self.p,mu_eV=0.)
        # finite numerical cancellation with symmetric valence/conduction
        self.assertLess(abs(eq20_physical(p)),abs(eq20_physical(self.p))*1e-6)

if __name__=='__main__':unittest.main()
