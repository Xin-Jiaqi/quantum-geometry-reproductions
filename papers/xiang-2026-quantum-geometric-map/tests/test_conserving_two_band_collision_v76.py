"""v76: true energy-shell two-band, CP dissipator, B-collision contact and two-color controls."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.conserving_two_band_collision_v76 import (
    Spec,CoherentSpec,shell,band_population,two_band_populations,
    population_weighted_spectral,bx_fixed_radius_nonelastic,coherent_data,
    apply_generator,vertex,_resolvent_spectral,_resolvent_dense,coherent_two_color,
    coherent_frame_covariance,coherent_static_spin,decoherence_nonuniqueness,
    equilibrium_B_collision_contact,mass_zeeman_collision_vertex,time_evolution,
    orbital_dissipator_mass_field,fermi,SX,SY,SZ)

class V76ConservingTwoBand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(n=16)
        cls.p=two_band_populations(cls.s)
        cls.c=CoherentSpec(n=16)
        cls.dc=coherent_data(cls.c)
        cls.ev=coherent_two_color(cls.c)
        cls.ind=coherent_two_color(cls.c,method='dense')
        cls.no=decoherence_nonuniqueness(cls.c)
        cls.Bz=mass_zeeman_collision_vertex(cls.c,delta=1e-5)

    def test_01_two_bands_present(self):self.assertEqual(set(self.p['bands']),{-1,1})
    def test_02_two_energies_each(self):
        for b in (-1,1):self.assertEqual(len(self.p['bands'][b]['entries']),2)
    def test_03_exact_both_band_isoenergy(self):
        for E in self.s.energies:
            for b in (-1,1):
                d=shell(self.s,E,b,B=.13)
                self.assertLess(np.max(abs(d['band_energy']-b*E)),2e-14)
    def test_04_Bx_rigid_translation_of_lab_ky(self):
        for b in (-1,1):
            a=shell(self.s,.3,b,0)
            c=shell(self.s,.3,b,.17)
            self.assertLess(np.max(abs(a['ky']-c['ky']-.17)),2e-14)
            np.testing.assert_allclose(a['kx'],c['kx'],atol=1e-14)
    def test_05_true_shell_DOS_weight(self):
        d=shell(self.s,.3,1)
        c=np.cos(d['theta']);target=(1+self.s.beta*c)**-2
        np.testing.assert_allclose(d['weights'],target/target.sum(),atol=1e-14)
    def test_06_band_reversal_DOS_not_same_at_each_angle(self):
        a=shell(self.s,.3,1);b=shell(self.s,.3,-1)
        self.assertGreater(np.linalg.norm(a['weights']-b['weights']),.04)
    def test_07_positive_transition_rates(self):
        for b in (-1,1):self.assertGreaterEqual(np.min(shell(self.s,.3,b)['transition']),0)
    def test_08_elastic_energy_spread_zero(self):
        for b in (-1,1):self.assertLess(np.ptp(shell(self.s,.6,b,.03)['band_energy']),2e-14)
    def test_09_weighted_detailed_balance(self):
        for b in (-1,1):
            d=shell(self.s,.3,b);w=d['weights'];W=d['transition']
            self.assertLess(np.max(abs(w[:,None]*W-w[None,:]*W.T)),2e-15)
    def test_10_charge_number_conservation_each_band(self):
        for b in (-1,1):
            d=shell(self.s,.3,b)
            self.assertLess(np.linalg.norm(d['weights']@d['L']),1e-14)
    def test_11_L_unity_zero_mode(self):
        for b in (-1,1):self.assertLess(np.linalg.norm(shell(self.s,.3,b)['L']@np.ones(self.s.n)),1e-14)
    def test_12_unique_positive_population_gap(self):
        for b in (-1,1):
            d=shell(self.s,.3,b);w=np.sqrt(d['weights'])
            ev=np.linalg.eigvalsh(w[:,None]*d['L']/w[None,:])
            self.assertLess(abs(ev[0]),1e-12)
            self.assertGreater(ev[1],.015)
    def test_13_single_band_population_source_conserved(self):
        for b in (-1,1):
            a=band_population(self.s,b)
            for e in a['entries']:
                self.assertLess(abs(e['source_mean']),1e-14)
                self.assertLess(abs(e['number']),1e-11)
    def test_14_direct_population_solve_residual(self):
        for b in (-1,1):
            for e in band_population(self.s,b)['entries']:
                self.assertLess(e['linear_residual'],1e-12)
    def test_15_independent_band_spectral(self):
        for b in (-1,1):
            self.assertLess(abs(band_population(self.s,b)['total_current']-
                population_weighted_spectral(self.s,b)),1e-10)
    def test_16_independent_ac_spectral(self):
        for b in (-1,1):
            self.assertLess(abs(band_population(self.s,b,.11)['total_current']-
                population_weighted_spectral(self.s,b,.11)),1e-10)
    def test_17_Bx_transport_rigid_translation(self):
        for B in (-.21,-.04,.08,.32):
            self.assertLess(abs(two_band_populations(self.s,B=B)['current']-self.p['current']),1e-12)
    def test_18_total_is_sum_bands(self):
        self.assertAlmostEqual(self.p['current'].real,10.962386926060933,places=10)
        self.assertLess(abs(self.p['current']-sum(b['total_current'] for b in self.p['bands'].values())),1e-13)
    def test_19_valence_nonzero_and_not_equal_conduction(self):
        self.assertGreater(abs(self.p['bands'][-1]['total_current']),1.)
        self.assertGreater(abs(self.p['bands'][1]['total_current']-self.p['bands'][-1]['total_current']),1.)
    def test_20_fixed_radius_not_elastic_for_Bx(self):
        a=bx_fixed_radius_nonelastic(self.s)
        self.assertTrue(a['fixed_radius_not_elastic'])
        self.assertGreater(a['conduction_spread'],.15)
    def test_21_untilted_fixed_radius_Bx_still_not_elastic(self):
        self.assertGreater(bx_fixed_radius_nonelastic(replace(self.s,beta=0.))['conduction_spread'],.1)
    def test_22_numerical_B_derivative_population_zero(self):
        h=.003
        v=lambda B:two_band_populations(self.s,B=B)['current']
        self.assertLess(abs((v(h)-v(-h))/(2*h)),1e-10)
    def test_23_coherence_h_band_energy_constant(self):
        self.assertLess(np.max(abs(self.dc['e'][:,0]+np.hypot(self.c.energy,self.c.mass))),1e-13)
        self.assertLess(np.ptp(self.dc['e'][:,1]),1e-13)
    def test_24_same_band_massive_spin_overlap(self):
        m=self.c.mass**2;r=self.c.energy**2
        t=self.dc['theta'];ratio=1+(r*np.cos(t[:,None]-t[None,:])+m)/(r+m)
        W=self.c.gamma/self.c.n*ratio;np.fill_diagonal(W,0)
        np.testing.assert_allclose(W,self.dc['W'],atol=1e-13)
    def test_25_coherence_GKSL_trace(self):
        rng=np.random.default_rng(81)
        X=rng.normal(size=(self.c.n,2,2))+1j*rng.normal(size=(self.c.n,2,2))
        self.assertLess(abs(np.trace(apply_generator(self.c,X),axis1=-2,axis2=-1).sum()),1e-13)
    def test_26_unperturbed_equilibrium_stationary(self):
        self.assertLess(np.linalg.norm(apply_generator(self.c,self.dc['rho0'])),1e-13)
    def test_27_extra_dephasing_leaves_populations(self):
        d=coherent_data(replace(self.c,pure_phi=.2))
        self.assertLess(np.linalg.norm(self.dc['L']-d['L']),1e-15)
        self.assertLess(np.linalg.norm(self.dc['W']-d['W']),1e-15)
    def test_28_generator_offdiagonal_dephasing_effect(self):
        X=np.zeros((self.c.n,2,2),complex);X[:,0,1]=.2
        a=apply_generator(self.c,X)
        b=apply_generator(replace(self.c,pure_phi=.2),X)
        self.assertLess(np.linalg.norm((b-a)[:,0,1]+.2*X[:,0,1]),1e-13)
    def test_29_positive_complete_GKSL_time_state(self):
        for t in (0,1,10,40):
            rho=time_evolution(self.c,t=t)
            self.assertGreater(np.min(np.linalg.eigvalsh(rho)),-1e-11)
            self.assertLess(abs(np.trace(rho,axis1=-2,axis2=-1).sum()-self.c.n*.93),1e-12)
    def test_30_GKSL_time_hermitian(self):
        rho=time_evolution(self.c,t=15)
        self.assertLess(np.linalg.norm(rho-rho.conj().swapaxes(-1,-2)),1e-13)
    def test_31_interband_source_from_commutator(self):
        source=vertex(self.dc['spin'],self.dc['rho0'])
        self.assertLess(np.max(abs(source[:,0,0])),1e-14)
        self.assertGreater(np.linalg.norm(source[:,0,1]),.2)
    def test_32_population_source_zero_for_coherence_only(self):
        source=vertex(self.dc['ry'],self.dc['rho0'])
        self.assertLess(np.max(abs(source[:,0,0])),1e-14)
    def test_33_dense_vs_spectral_two_colors(self):
        for name in ('b_then_e','e_then_b','total','first_e','first_b','local_harmonic_1'):
            self.assertLess(abs(self.ev[name]-self.ind[name]),2e-12,name)
    def test_34_both_mixed_orders_present_locally(self):
        self.assertGreater(self.ev['local_mixed_norm'],5)
        self.assertGreater(abs(self.ev['local_harmonic_1']),.3)
    def test_35_angular_integrated_EB_selection_zero(self):
        self.assertLess(abs(self.ev['total']),2e-12)
        self.assertIn('NOT_AUTHOR_LENGTH_GAUGE',self.ev['status'])
    def test_36_trace_two_color_Ward(self):
        for name in ('trace_first_e','trace_first_b','trace_second'):
            self.assertLess(abs(self.ev[name]),1e-12)
    def test_37_two_color_ordering_exact(self):
        self.assertLess(self.ev['two_color_order_residual'],1e-12)
    def test_38_nontrivial_momentum_U1_frame(self):
        g=coherent_frame_covariance(self.c,.47)
        for name in g:self.assertLess(abs(g[name]-self.ev[name]),2e-12,name)
    def test_39_second_U1_frame(self):
        g=coherent_frame_covariance(self.c,-.91)
        self.assertLess(abs(g['local_harmonic_1']-self.ev['local_harmonic_1']),1e-12)
    def test_40_same_transport_collision_not_unique_coherence(self):
        self.assertLess(self.no['population_L_error'],1e-15)
        self.assertGreater(abs(self.no['spin_difference']),.05)
    def test_41_same_tau_distinct_mixed_local_response(self):
        self.assertGreater(abs(self.no['local_harmonic_difference']),.15)
    def test_42_Bx_equilibrium_missing_collision(self):
        v=equilibrium_B_collision_contact(self.c,B=1e-4)
        self.assertLess(v['analytic_coherent_identity_residual'],1e-12)
        self.assertGreater(v['frozen_B_stationarity_defect'],.2)
    def test_43_Bx_necessary_contact_trace_conserved(self):
        v=equilibrium_B_collision_contact(self.c,B=1e-4)
        self.assertLess(abs(v['required_contact_total_trace']),1e-12)
        self.assertLess(v['finite_B_derivative_residual'],2e-7)
    def test_44_true_Bz_Lindblad_collision_eq(self):
        self.assertLess(self.Bz['field_equilibrium_collision_defect_at_zero'],1e-12)
        self.assertLess(self.Bz['equilibrium_collision_ward_residual'],5e-9)
    def test_45_Bz_magnetic_vertex_essential(self):
        self.assertGreater(self.Bz['missing_collision_vertex_error'],.07)
        self.assertGreater(self.Bz['collision_vertex_drive_norm'],.3)
    def test_46_Bz_full_response_five_point(self):
        self.assertLess(self.Bz['full_difference_error'],2e-9)
    def test_47_Bz_field_input_distinct_source_Bx(self):
        self.assertIn('NOT_ORIGINAL_BX_PHE',self.Bz['source'])
        self.assertIn('NOT_CONSTRUCTED_FULL_B_DEPENDENT_BATH',equilibrium_B_collision_contact(self.c)['status'])
    def test_48_Bz_limit_step_halves_error(self):
        a=mass_zeeman_collision_vertex(self.c,delta=5e-5)
        self.assertLess(self.Bz['full_difference_error'],a['full_difference_error']*.2)
    def test_49_spec_fail_closed(self):
        for kw in ({'n':13},{'gamma':0},{'temperature':0},{'beta':1.},{'energies':(-.1,)},{'eta_b':0}):
            with self.assertRaises(ValueError):Spec(**kw)
    def test_50_coherent_spec_fail_closed(self):
        for kw in ({'n':9},{'pure_phi':-.1},{'gamma':0},{'eta_e':0},{'energy':0}):
            with self.assertRaises(ValueError):CoherentSpec(**kw)
    def test_51_unknown_method_and_negative_times(self):
        with self.assertRaises(ValueError):coherent_two_color(self.c,'wrong')
        with self.assertRaises(ValueError):time_evolution(self.c,t=-2)
        with self.assertRaises(ValueError):mass_zeeman_collision_vertex(self.c,delta=-1)
    def test_52_population_AC_preserves_B_covariance(self):
        for b in (-1,1):
            self.assertLess(abs(band_population(self.s,b,omega=.1,B=.2)['total_current']-
                band_population(self.s,b,omega=.1,B=0)['total_current']),1e-12)

if __name__=='__main__':unittest.main()
