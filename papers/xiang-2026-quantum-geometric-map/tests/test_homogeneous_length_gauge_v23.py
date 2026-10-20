"""v23: homogeneous length-gauge model gates, NOT experimental PHE validation."""
import unittest
import numpy as np
from dataclasses import replace
from qgmap.homogeneous_length_gauge_v23 import (UniformSpec,kmesh,h0_and_vel,inputs,fft_dy,cov_dy,
    orb_to_band,fermi,source_and_channels,algebraic_channels,leading_fermi_pole,eq20_direct_band_current,
    prepared_slope,population_preparation_contribution,comm)

class TestHomogeneousV23(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=UniformSpec(n=64,omega_e=.023,omega_b=.015,eta_e=.07,eta_b=.04)
        cls.fft=source_and_channels(cls.p)
        cls.alg=algebraic_channels(cls.p)
    def test_01_periodic_fft_derivative(self):
        p=self.p;k=kmesh(p.n); a=np.cos(3*k[...,1]); self.assertLess(np.max(np.abs(fft_dy(a)+3*np.sin(3*k[...,1]))),1e-12)
    def test_02_band_gap_not_degenerate(self):
        self.assertGreater(self.fft['gap_min'],.5)
    def test_03_e_then_b_independent_derivative(self):
        self.assertLess(abs(self.fft['e_then_b']-self.alg['e_then_b']),3e-8)
    def test_04_b_then_e_independent_derivative(self):
        self.assertLess(abs(self.fft['b_then_e']-self.alg['b_then_e']),3e-8)
    def test_05_full_two_orderings(self):
        self.assertLess(abs(self.fft['total']-self.fft['b_then_e']-self.fft['e_then_b']),3e-14)
        self.assertLess(abs(self.fft['total']-self.alg['total']),3e-8)
    def test_06_density_trace(self):
        self.assertLess(abs(np.mean(np.trace(self.fft['rho_eb'],axis1=-2,axis2=-1))),1e-12)
    def test_07_zero_equilibrium_magnetic_population(self):
        h,_,_,_,_=inputs(self.p)
        _,u=np.linalg.eigh(h)
        db=orb_to_band(u,self.fft['rho_b'])
        self.assertLess(np.max(abs(np.diagonal(db,axis1=-2,axis2=-1))),1e-13)
    def test_08_eq20_pole_independent(self):
        q=UniformSpec(n=64)
        self.assertLess(abs(leading_fermi_pole(q)-eq20_direct_band_current(q)),1e-13)
    def test_09_adiabatic_pole(self):
        base=UniformSpec(n=64)
        goal=eq20_direct_band_current(base)
        a=replace(base,eta_e=.02,eta_b=.02)
        b=replace(base,eta_e=.005,eta_b=.005)
        erra=abs(.02*algebraic_channels(a)['total']-goal)
        errb=abs(.005*algebraic_channels(b)['total']-goal)
        self.assertLess(errb,4e-7)
        self.assertLess(errb,erra*.15)
    def test_10_fixed_real_frequency_no_simple_tau_pole(self):
        base=UniformSpec(n=64,omega_e=.25,omega_b=.1)
        a=replace(base,eta_e=.04,eta_b=.04)
        b=replace(base,eta_e=.01,eta_b=.01)
        x=.04*algebraic_channels(a)['total'];y=.01*algebraic_channels(b)['total']
        self.assertLess(abs(y),.5*abs(x))
    def test_11_static_prepared_population_identity(self):
        p=UniformSpec(n=64,omega_e=.12,eta_e=.09,omega_b=0.)
        full=prepared_slope(p,b=.0001)
        frozen=prepared_slope(p,b=.0001,frozen=True)
        pop=population_preparation_contribution(p)
        self.assertLess(abs(full-frozen-pop),2e-7)
        self.assertGreater(abs(pop),.005)
    def test_12_dynamic_magnetic_switch_matches_frozen_limit(self):
        p=UniformSpec(n=64,omega_e=.12,eta_e=.09,omega_b=0.,eta_b=.0001)
        frozen=prepared_slope(p,b=.0001,frozen=True)
        dyn=algebraic_channels(p)['total']
        self.assertLess(abs(dyn-frozen),2e-6)
    def test_13_gauge_covariant_resolvent(self):
        other=replace(self.p,gauge_amplitude=.31)
        al=algebraic_channels(other)
        for k in ('b_then_e','e_then_b','total'):
            self.assertLess(abs(al[k]-self.alg[k]),1e-12)
    def test_14_gauge_covariant_fft_discretization(self):
        other=replace(self.p,gauge_amplitude=.31)
        z=source_and_channels(other)
        self.assertLess(abs(z['total']-self.fft['total']),5e-9)
    def test_15_covariant_derivative_of_zeeman_zero(self):
        p=replace(self.p,gauge_amplitude=.31)
        h,vx,vy,z,a=inputs(p)
        self.assertLess(np.max(abs(cov_dy(z,a))),1e-12)
    def test_16_input_parameters_fail_closed(self):
        with self.assertRaises(ValueError):UniformSpec(n=7)
        with self.assertRaises(ValueError):UniformSpec(temp=0)
        with self.assertRaises(ValueError):UniformSpec(eta_e=0)
    def test_17_band_derivatives_correct(self):
        p=self.p; h,vx,vy=h0_and_vel(p)
        self.assertLess(np.max(abs(fft_dy(h)-vy)),2e-12)
    def test_18_k_mesh_convergence(self):
        p=replace(self.p,n=48)
        x=algebraic_channels(p)['total']
        self.assertLess(abs(x-self.alg['total']),3e-6)

if __name__=='__main__':unittest.main()
