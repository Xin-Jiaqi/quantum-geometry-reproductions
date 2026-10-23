"""v74: physical source Eq15 pole-residue, two-order cancellation and limit tests."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.selected_electric_pole_v74 import (
    Spec, channels, paths, finite_frequency_negative_control, literal_eq20,
    printed_eq31, _data, _resolv, _comm, SX, _integrate)

class V74SelectedPole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(nr=64,nphi=88,eta_e=.003,eta_b=.001)
        cls.v=channels(cls.s)
        cls.literal=literal_eq20(cls.s)
        cls.curves=paths(cls.s)
        cls.frequency=finite_frequency_negative_control(cls.s)
        cls.sm=channels(replace(cls.s,helicity='sm'))
        cls.xy=channels(cls.s)
        cls.yx=channels(replace(cls.s,direction='yx'))

    def test_01_full_is_sum_of_both_insertions(self):
        v=self.v;self.assertLess(abs(v['full']-v['b_then_e']-v['e_then_b']),1e-13)

    def test_02_e_then_b_is_exact_intraband_plus_other(self):
        self.assertLess(self.v['eb_decomposition_error'],1e-12)

    def test_03_intraband_E_pole_is_nonzero(self):
        self.assertGreater(abs(self.v['eta_e_scaled_intraband_e_pole']),.002)

    def test_04_B_then_E_is_not_negligible(self):
        self.assertGreater(abs(self.v['eta_e_scaled_b_then_e']),.003)

    def test_05_E_then_B_is_not_negligible(self):
        self.assertGreater(abs(self.v['eta_e_scaled_e_then_b']),.006)

    def test_06_two_orders_cancel_large_contributions(self):
        self.assertLess(abs(self.v['eta_e_scaled_full']),.45*abs(self.v['eta_e_scaled_e_then_b']))

    def test_07_full_residue_near_source_literal_eq20(self):
        self.assertLess(abs(self.v['eta_e_scaled_full']-self.literal),5e-7)

    def test_08_selective_Ediag_residue_near_source_literal(self):
        self.assertLess(abs(self.v['eta_e_scaled_intraband_e_pole']-self.literal),2e-7)

    def test_09_selected_pole_residue_not_full_unregulated_current(self):
        self.assertGreater(abs(self.v['full']-self.literal),.2)

    def test_10_printed_eq31_still_opposite(self):
        self.assertLess(abs(printed_eq31(self.s)+self.literal),1e-14)
        self.assertGreater(abs(self.v['eta_e_scaled_full']-printed_eq31(self.s)),.0049)

    def test_11_source_literal_sign_not_fitted(self):
        self.assertLess(self.literal.real,-.002)
        self.assertGreater(printed_eq31(self.s).real,.002)

    def test_12_regulator_sum_is_preserved(self):
        self.assertLess(abs(self.v['output_regulator']-.004),1e-14)

    def test_13_exact_analytic_fermi_derivative_frechet(self):
        d=_data(self.s)
        gap=d['e'][:,:,None]-d['e'][:,None,:]
        norm=np.linalg.norm(d['De']-_resolv(d,d['De'],1.+0j)-1j*(_comm(d['h'],_resolv(d,d['De'],1.+0j))))
        self.assertLess(norm,1e-10)

    def test_14_Ddiag_is_exact_intraband(self):
        d=_data(self.s)
        val=d['ud']@d['Ddiag']@d['u']
        diagfree=val.copy()
        idx=np.arange(2);diagfree[...,idx,idx]=0
        self.assertLess(np.max(abs(diagfree)),2e-13)

    def test_15_Doff_has_zero_intraband_component(self):
        d=_data(self.s)
        val=d['ud']@d['Doff']@d['u']
        self.assertLess(np.max(abs(np.diagonal(val,axis1=-2,axis2=-1))),2e-13)

    def test_16_analytic_frechet_matches_rho0_finite_y_difference(self):
        d=_data(self.s)
        x,y=d['x'],d['y'];beta=self.s.beta;sg=1
        def matrix(yy):
            h=beta*x[:,None,None]*np.eye(2)+sg*(x[:,None,None]*np.array([[0,-1j],[1j,0]])-yy[:,None,None]*SX)
            e,u=np.linalg.eigh(h)
            f=1/(1+np.exp(np.clip((e-self.s.mu)/self.s.temperature,-600,600)))
            return (u*f[:,None,:])@np.swapaxes(u.conj(),-1,-2)
        step=2e-5
        fd=(matrix(y+step)-matrix(y-step))/(2*step)
        self.assertLess(np.max(abs(fd-d['De'])),3e-7)

    def test_17_electric_offdiagonal_discard_breaks_full_current(self):
        self.assertGreater(abs(self.v['full']-self.v['electric_diag_pole']),3e-5)

    def test_18_unselected_terms_cancel_not_individually_zero(self):
        self.assertGreater(abs(self.v['b_then_e']),1.)
        self.assertGreater(abs(self.v['eb_non_electric_diag']),1.)
        self.assertLess(abs(self.v['b_then_e']+self.v['eb_non_electric_diag']),.002)

    def test_19_balanced_eta_path_selected_residue_decreases_error(self):
        rows=self.curves['etaB_equals_etaE']
        errors=[abs(x['scaled_full']-self.literal) for x in rows]
        self.assertLess(errors[-1],errors[0]/20)

    def test_20_balanced_eta_path_final_eq20_error(self):
        self.assertLess(abs(self.curves['etaB_equals_etaE'][-1]['scaled_full']-self.literal),2e-7)

    def test_21_etaB_tenth_path_matches_source_residue(self):
        self.assertLess(abs(self.curves['etaB_tenth_etaE'][-1]['scaled_full']-self.literal),2e-7)

    def test_22_etaB_ten_path_still_approaches(self):
        self.assertLess(abs(self.curves['etaB_ten_etaE'][-1]['scaled_full']-self.literal),2e-6)

    def test_23_individual_orderings_path_dependent(self):
        a=self.curves['etaB_tenth_etaE'][0]['scaled_B_then_E']
        b=self.curves['etaB_ten_etaE'][0]['scaled_B_then_E']
        self.assertGreater(abs(a-b),.001)

    def test_24_full_combined_residue_less_path_dependent(self):
        vals=[r[-1]['scaled_full'] for r in self.curves.values()]
        self.assertLess(max(abs(x-vals[0]) for x in vals),3e-6)

    def test_25_finite_omega_E_scaled_full_vanishes(self):
        self.assertLess(abs(self.frequency[-1][1]),.2*abs(self.frequency[0][1]))

    def test_26_finite_omega_E_selective_pole_residue_vanishes(self):
        self.assertLess(abs(self.frequency[-1][2]),.2*abs(self.frequency[0][2]))

    def test_27_finite_real_frequency_no_wrong_tau_label(self):
        self.assertGreater(abs(self.frequency[-1][1]),1e-6)
        self.assertLess(abs(self.frequency[-1][1]),.00002)

    def test_28_main_and_sm_helicity_opposite(self):
        self.assertLess(abs(self.sm['eta_e_scaled_full']+self.v['eta_e_scaled_full']),1e-12)

    def test_29_main_and_sm_literal_opposite(self):
        self.assertLess(abs(literal_eq20(self.s)+literal_eq20(replace(self.s,helicity='sm'))),1e-12)

    def test_30_both_tensor_directions_finite(self):
        self.assertTrue(np.isfinite([self.xy['full'].real,self.yx['full'].real]).all())

    def test_31_tensor_yx_nonzero(self):
        self.assertGreater(abs(self.yx['eta_e_scaled_full']),.0001)

    def test_32_electric_source_pole_is_selected_not_global_gamma(self):
        self.assertIn('SELECTED_E_POLE',self.v['status'])
        self.assertIn('NOT_FULL_DC_PHE',self.v['status'])

    def test_33_magnetic_omega_negative_real_amp_response(self):
        a=channels(replace(self.s,omega_e=.1,omega_b=-.03))['full']
        b=channels(replace(self.s,omega_e=-.1,omega_b=.03))['full']
        self.assertLess(abs(a-b.conjugate()),1e-9)

    def test_34_etaB_fixed_large_does_not_establish_zero_B_pole_limit(self):
        a=channels(replace(self.s,eta_e=.001,eta_b=.05))
        self.assertGreater(abs(a['eta_e_scaled_full']-self.literal),1e-6)

    def test_35_mu_changes_literal_pole(self):
        a=literal_eq20(replace(self.s,mu=-.3))
        self.assertLess(abs(a+self.literal),2e-12)
        self.assertGreater(abs(a-self.literal),.004)

    def test_36_negative_frequency_is_not_reversed_eta(self):
        a=channels(replace(self.s,omega_e=.08,omega_b=.02))
        b=channels(replace(self.s,omega_e=-.08,omega_b=-.02))
        self.assertLess(abs(a['full']-b['full'].conjugate()),1e-9)

    def test_37_positive_regulator_fails_closed(self):
        for kw in ({'eta_e':0},{'eta_b':0},{'temperature':0}):
            with self.assertRaises(ValueError):Spec(**kw)

    def test_38_helicity_fails_closed(self):
        with self.assertRaises(ValueError):Spec(helicity='uncertain')

    def test_39_direction_fails_closed(self):
        with self.assertRaises(ValueError):Spec(direction='xx')

    def test_40_nonfinite_source_fails_closed(self):
        with self.assertRaises(ValueError):Spec(omega_e=np.inf)
        with self.assertRaises(ValueError):Spec(beta=np.nan)

    def test_41_quadrature_positive_and_no_degeneracy(self):
        d=_data(self.s)
        self.assertGreater(min(d['weight']),0)
        self.assertGreater(np.min(abs(d['e'][:,0]-d['e'][:,1])),1e-7)

    def test_42_source_eq20_is_real_physical_B_coefficient(self):
        self.assertLess(abs(self.literal.imag),1e-12)

    def test_43_eta_scaled_full_real_on_zero_freq(self):
        self.assertLess(abs(self.v['eta_e_scaled_full'].imag),1e-10)

    def test_44_exact_resolvent_generic_complex_source(self):
        d=_data(self.s);rng=np.random.default_rng(21)
        src=rng.normal(size=d['h'].shape)+1j*rng.normal(size=d['h'].shape)
        z=.023-.11j;r=_resolv(d,src,z)
        self.assertLess(np.max(abs(z*r+1j*_comm(d['h'],r)-src)),1e-13)

    def test_45_source_analytic_no_fitted_sign_or_factor(self):
        self.assertLess(abs(self.v['eta_e_scaled_full']-self.literal),1e-6)
        self.assertGreater(abs(self.v['eta_e_scaled_full']+self.literal),.004)

if __name__=='__main__':unittest.main()
