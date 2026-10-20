"""v25 reproducibility tests. Paper-figure discrepancies are audited, not 'fixed'."""
import math
import unittest
import numpy as np
from qgmap.si_fig2_v25 import (FigureInputs, H, HBAR, E, MU_B, geometry, geometry_printed,
                               independent_angular_geometry, tanh_population, evaluate,
                               independent_volt_from_ratio)

class TestFigureSI(unittest.TestCase):
    def setUp(self): self.p=FigureInputs()
    def test_stable_closed_form_matches_printed(self):
        for beta in (.01,.0625,.125,.25,.5,.8):
            self.assertAlmostEqual(geometry(beta),geometry_printed(beta),delta=2e-12)
    def test_direct_angular_quadrature(self):
        for beta in (.05,.125,.3,.7):
            self.assertLess(abs(geometry(beta)-independent_angular_geometry(beta,2048)),1e-13)
    def test_zero_t_figure_geometry(self):
        self.assertAlmostEqual(geometry(self.p.beta),-.00250641590669,delta=1e-12)
    def test_ordinate_si(self):
        self.assertAlmostEqual(evaluate(self.p)['minus_1e4_sigma_over_e2_hbar'],11.02083611497046,delta=1e-9)
    def test_eq11_voltage(self):
        self.assertAlmostEqual(evaluate(self.p)['open_circuit_voltage_microV'],268.26200776,delta=1e-7)
    def test_2pi_diagnostic_ordinate(self):
        q=evaluate(self.p,hypothesis='single_2pi_diagnostic')
        self.assertAlmostEqual(q['minus_1e4_sigma_over_e2_hbar'],1.75402054470,delta=1e-10)
    def test_2pi_diagnostic_voltage(self):
        self.assertAlmostEqual(evaluate(self.p,hypothesis='single_2pi_diagnostic')['open_circuit_voltage_microV'],42.695224579,delta=1e-7)
    def test_ratio_exactly_2pi(self):
        x=evaluate(self.p)['sigma_sheet_S'];y=evaluate(self.p,hypothesis='single_2pi_diagnostic')['sigma_sheet_S']
        self.assertAlmostEqual(x/y,2*math.pi,delta=1e-13)
    def test_independent_conductance_quantum_path(self):
        for useh,branch in ((False,'literal_hbar'),(True,'single_2pi_diagnostic')):
            a=evaluate(self.p,hypothesis=branch)['open_circuit_voltage_V']
            b=independent_volt_from_ratio(self.p,h_use=useh)
            self.assertLess(abs(a-b),5e-15)
    def test_g_tau_field_linearity(self):
        from dataclasses import replace
        p=self.p
        s=evaluate(p)['sigma_sheet_S']
        for kwargs in ({'g':p.g*3},{'tau_s':p.tau_s*3},{'B_tesla':p.B_tesla*3}):
            self.assertAlmostEqual(evaluate(replace(p,**kwargs))['sigma_sheet_S']/s,3.,delta=1e-13)
    def test_mu_parity_and_zero(self):
        for temp in (0.,20.,40.):
            a=evaluate(self.p,20.,temp)['sigma_sheet_S']
            b=evaluate(self.p,-20.,temp)['sigma_sheet_S']
            z=evaluate(self.p,0.,temp)['sigma_sheet_S']
            self.assertAlmostEqual(a,-b,delta=1e-16)
            self.assertEqual(z,0.)
    def test_step_thermal_broadening(self):
        for mu in (5.,10.,20.,40.):
            vals=[abs(evaluate(self.p,mu,T)['sigma_sheet_S']) for T in (0,20,40)]
            self.assertGreater(vals[0],vals[1]);self.assertGreater(vals[1],vals[2])
    def test_sheet_conductivity_and_voltage_dimensions(self):
        # A transverse conductivity in S gives V = S*Ohm*(V/m)*m.
        a=evaluate(self.p)
        self.assertAlmostEqual(a['open_circuit_voltage_V'],
                               -a['sigma_per_B_S_per_T']*self.p.B_tesla*self.p.resistance_ohms*self.p.E_vm*self.p.transverse_width_m)
    def test_ambiguous_unit_hypothesis_is_explicit(self):
        with self.assertRaises(ValueError):evaluate(self.p,hypothesis='undocumented_fix')
        with self.assertRaises(ValueError):geometry(1.1)
    def test_density_dimensionless_comparison_fig(self):
        # Plateau from original figure is visually ~1.75 (not digitized data).
        fig_estimate=1.75
        predicted=evaluate(self.p)['minus_1e4_sigma_over_e2_hbar']
        alt=evaluate(self.p,hypothesis='single_2pi_diagnostic')['minus_1e4_sigma_over_e2_hbar']
        self.assertGreater(abs(predicted-fig_estimate),9.)
        self.assertLess(abs(alt-fig_estimate),0.04)

if __name__=='__main__':unittest.main()
