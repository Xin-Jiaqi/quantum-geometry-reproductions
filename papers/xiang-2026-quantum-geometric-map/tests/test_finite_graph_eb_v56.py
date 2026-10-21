"""v56 independent finite E/B response tests, NOT a bulk orbital MNHE claim."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.finite_graph_eb_v56 import (
    Spec,fermi,equilibrium_density,frechet_fermi,polarization,current,resolvent,
    state,independent_sigma_spectral,mixed_eb,finite_field_derivative,equilibrium_current)
from qgmap.embedding_magnetic_v55 import (peierls_field,field_gauge_unitary)

class FiniteGraphEBV56(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=np.array([[0,0],[1.25,0],[1.25,.8],[0,.8]],float)
        cls.h=np.array([[.2,1.,.3,.8],[1.,-.1,1.2,.1],[.3,1.2,.05,.9],[.8,.1,.9,.3]],complex)
        cls.s=Spec()
        cls.ref=mixed_eb(cls.h,cls.p,cls.s)
        cls.land=mixed_eb(cls.h,cls.p,cls.s,gauge='landau_y')

    def test_01_original_peierls_field_and_signed_charge_propagation(self):
        p=peierls_field(self.h,self.p,self.s.field,self.s.q_over_hbar)
        z=state(self.h,self.p,self.s)
        np.testing.assert_allclose(z['H'],p.hamiltonian,atol=1e-14)
        np.testing.assert_allclose(z['M'],p.magnetic_vertex,atol=1e-14)

    def test_02_static_electric_polarization_commutator_current(self):
        z=state(self.h,self.p,self.s)
        np.testing.assert_allclose(z['J'],1j*(z['H']@z['P_current']-z['P_current']@z['H']),atol=1e-14)
        self.assertAlmostEqual(np.trace(z['J']),0,delta=1e-14)
        np.testing.assert_allclose(z['J'],z['J'].conj().T,atol=1e-14)

    def test_03_spectral_sum_independent_of_liouville_linear_solver(self):
        z=state(self.h,self.p,self.s)
        x=independent_sigma_spectral(z['H'],z['P_current'],z['P_drive'],self.s.mu,self.s.temperature,self.s.gamma,self.s.omega)
        self.assertAlmostEqual(z['sigma'],x,delta=1e-12)
        self.assertGreater(abs(z['sigma']),.005)

    def test_04_liouville_equation_matrix_residual(self):
        z=state(self.h,self.p,self.s)
        lhs=z['z']*z['rho_E']+1j*(z['H']@z['rho_E']-z['rho_E']@z['H'])
        np.testing.assert_allclose(lhs,z['source_E'],atol=1e-13)

    def test_05_exact_two_magnetic_gauges_electric_sigma(self):
        self.assertAlmostEqual(self.ref['response_sigma'],self.land['response_sigma'],delta=3e-13)

    def test_06_analytic_complete_mixed_B_vs_independent_five_point_B_difference(self):
        for gauge in ('symmetric','landau_y'):
            d=mixed_eb(self.h,self.p,self.s,gauge)
            finite=finite_field_derivative(self.h,self.p,self.s,gauge,step=1e-4)
            self.assertAlmostEqual(d['mixed_dsigma_dB'],finite,delta=2e-10)

    def test_07_mixed_B_response_gauge_invariance_after_full_sum(self):
        self.assertAlmostEqual(self.ref['mixed_dsigma_dB'],self.land['mixed_dsigma_dB'],delta=1e-12)
        self.assertGreater(abs(self.ref['mixed_dsigma_dB']),.03)

    def test_08_components_individually_not_invariant_under_B_dependent_gauge(self):
        a=self.ref['contributions'];b=self.land['contributions']
        for term in a:
            self.assertGreater(abs(a[term]-b[term]),.005)
        self.assertAlmostEqual(sum(a.values()),sum(b.values()),delta=1e-12)

    def test_09_mandatory_current_contact_term_nonzero_and_omission_fails_B_derivative(self):
        a=self.ref['contributions']
        without=a['equilibrium_fermi_derivative']+a['magnetic_liouvillian']
        self.assertGreater(abs(without-self.ref['mixed_dsigma_dB']),.02)
        self.assertGreater(abs(without-finite_field_derivative(self.h,self.p,self.s)),.02)

    def test_10_complete_equilibrium_rho_B_Frechet_vs_independent_field_FD(self):
        for gauge in ('symmetric','landau_y'):
            z=state(self.h,self.p,self.s,gauge)
            e=1e-4
            plus=peierls_field(self.h,self.p,self.s.field+e,self.s.q_over_hbar,gauge=gauge).hamiltonian
            minus=peierls_field(self.h,self.p,self.s.field-e,self.s.q_over_hbar,gauge=gauge).hamiltonian
            ref=(equilibrium_density(plus,self.s.mu,self.s.temperature)-equilibrium_density(minus,self.s.mu,self.s.temperature))/(2*e)
            np.testing.assert_allclose(z['rho_B'],ref,atol=5e-10)

    def test_11_Frechet_degenerate_spectrum_gauge_invariant_and_FD(self):
        h=np.diag([0.,0.,1.,2.]).astype(complex)
        m=np.array([[.1,.24,.2j,0],[.24,-.2,0,-.2j],[-.2j,0,.23,.1],[0,.2j,.1,.3]])
        analytic=frechet_fermi(h,m,.2,.4)
        delta=2e-5
        finite=(equilibrium_density(h+delta*m,.2,.4)-equilibrium_density(h-delta*m,.2,.4))/(2*delta)
        np.testing.assert_allclose(analytic,finite,atol=1e-9)

    def test_12_frequency_reflection_is_complex_conjugate(self):
        a=mixed_eb(self.h,self.p,self.s)
        b=mixed_eb(self.h,self.p,replace(self.s,omega=-self.s.omega))
        self.assertAlmostEqual(a['response_sigma'],b['response_sigma'].conjugate(),delta=2e-13)
        self.assertAlmostEqual(a['mixed_dsigma_dB'],b['mixed_dsigma_dB'].conjugate(),delta=2e-13)

    def test_13_spinless_onsager_reciprocity_across_real_B_reversal(self):
        for B in (0.,.21,-.27):
            p=replace(self.s,field=B,probe_axis='x',drive_axis='y')
            q=replace(self.s,field=-B,probe_axis='y',drive_axis='x')
            self.assertAlmostEqual(state(self.h,self.p,p)['sigma'],state(self.h,self.p,q)['sigma'],delta=5e-13)
        at0=mixed_eb(self.h,self.p,replace(self.s,field=0.,probe_axis='x',drive_axis='y'))['mixed_dsigma_dB']
        reversed_index=mixed_eb(self.h,self.p,replace(self.s,field=0.,probe_axis='y',drive_axis='x'))['mixed_dsigma_dB']
        self.assertAlmostEqual(at0,-reversed_index,delta=5e-13)

    def test_14_charge_flip_with_B_flip_preserves_sigma_and_reverses_mixed_slope(self):
        a=self.s
        b=replace(a,q_over_hbar=-a.q_over_hbar,field=-a.field)
        sa=mixed_eb(self.h,self.p,a);sb=mixed_eb(self.h,self.p,b)
        self.assertAlmostEqual(sa['response_sigma'],sb['response_sigma'],delta=3e-13)
        self.assertAlmostEqual(sa['mixed_dsigma_dB'],-sb['mixed_dsigma_dB'],delta=3e-13)

    def test_15_closed_finite_equilibrium_dc_current_zero_even_at_static_E_B(self):
        for E in (0,.02,.2):
            for gauge in ('symmetric','landau_y'):
                z=equilibrium_current(self.h,self.p,self.s,gauge,E)
                self.assertLess(abs(z),3e-13)

    def test_16_zero_frequency_finite_closed_graph_gamma_to_zero_vanishes(self):
        a=abs(state(self.h,self.p,replace(self.s,omega=0,gamma=.001))['sigma'])
        b=abs(state(self.h,self.p,replace(self.s,omega=0,gamma=.0005))['sigma'])
        self.assertLess(b,a*.52)
        self.assertGreater(b,a*.49)

    def test_17_origin_shift_does_not_change_dipole_linear_response(self):
        for x,y in ((2.3,-1.7),(-.3,.4)):
            a=mixed_eb(self.h,self.p+np.array([x,y]),self.s)
            self.assertAlmostEqual(a['response_sigma'],self.ref['response_sigma'],delta=5e-13)
            self.assertAlmostEqual(a['mixed_dsigma_dB'],self.ref['mixed_dsigma_dB'],delta=5e-13)

    def test_18_fail_closed_bad_tau_T_charge_geometry_and_stencil(self):
        for vals in ({'temperature':0},{'gamma':0},{'q_over_hbar':0},{'probe_axis':'z'},{'field':float('nan')}):
            with self.assertRaises(ValueError):Spec(**vals)
        with self.assertRaises(ValueError):polarization(self.p,-1,'z')
        with self.assertRaises(ValueError):finite_field_derivative(self.h,self.p,self.s,step=0)
        with self.assertRaises(ValueError):finite_field_derivative(self.h,self.p,self.s,points=7)
        with self.assertRaises(ValueError):resolvent(np.diag([0,1]),np.eye(2),0j)

    def test_19_output_physics_status_is_conditional_not_whole_mnhe(self):
        self.assertIn('CONDITIONAL_FINITE_OPEN_GRAPH',self.ref['acceptance'])
        self.assertNotIn('VERIFIED_BULK',self.ref['acceptance'])
        self.assertIn('Not bulk MNHE',self.ref['limitation'])

if __name__=='__main__':unittest.main()
