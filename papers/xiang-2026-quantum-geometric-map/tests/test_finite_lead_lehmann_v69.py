"""v69 finite-Hamiltonian Lehmann oracle and ordered-limit science gates."""
import json
import pathlib
import unittest
from dataclasses import replace
import numpy as np
from qgmap.finite_lead_lehmann_v69 import (
    Spec, build_finite_system, lehmann_response, static_isothermal_response,
    equilibrium_current, static_finite_difference, finite_difference_B,
    finite_difference_static_mixed, lehmann_continuity_ward,
)
from qgmap.lead_keldysh_eb_v68 import Spec as KSpec, response as keldysh_response
from qgmap.dynamic_orbital_eb_v66 import comm

ROOT=pathlib.Path(__file__).resolve().parents[1]
FROZEN=json.loads((ROOT/'results/frozen/v69_finite_lead_lehmann.json').read_text())


class V69FiniteLeadLehmann(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(length=60,eta=.03)
        cls.frames={f:lehmann_response(cls.s,f) for f in ('symmetric','landau_y','orbital')}
        cls.medium=lehmann_response(Spec(length=250,eta=.02))
        cls.long=lehmann_response(Spec(length=350,eta=.02))
        cls.intermediate=lehmann_response(Spec(length=500,eta=.01))
        cls.production=lehmann_response(Spec(length=1000,eta=.005))
        cls.keldysh=keldysh_response(KSpec(nodes_per_interval=640))
        cls.long_B=finite_difference_B(Spec(length=350,eta=.02),step=.003,points=5)
        cls.static_spec=Spec(length=50,eta=.02,omega=0.)
        cls.static=static_isothermal_response(cls.static_spec)['sigma_device']
        cls.static_fd=static_finite_difference(cls.static_spec,step=1e-4)
        cls.ward=lehmann_continuity_ward(Spec(length=40,eta=.03))

    def test_01_full_hamiltonian_is_hermitian(self):
        b=build_finite_system(self.s)
        self.assertLess(np.linalg.norm(b['H']-b['H'].conj().T),1e-12)

    def test_02_full_hamiltonian_dimension(self):
        b=build_finite_system(Spec(length=31))
        self.assertEqual(b['H'].shape,(66,66))

    def test_03_exact_two_finite_leads(self):
        b=build_finite_system(Spec(length=31))
        self.assertEqual(b['length'],31)
        self.assertEqual(b['device_dimension'],4)

    def test_04_internal_lead_hopping(self):
        t=1.1;L=25;b=build_finite_system(Spec(length=L));H=b['H']
        self.assertAlmostEqual(H[4,5].real,t,14)
        self.assertAlmostEqual(H[4+L,5+L].real,t,14)

    def test_05_endpoint_contact_is_not_multiplied_by_t(self):
        b=build_finite_system(Spec(length=30));o=b['operators']
        self.assertLess(np.linalg.norm(b['H'][:4,4]-o['vectors'][0]),1e-14)
        self.assertGreater(np.linalg.norm(o['vectors'][0]),.1)

    def test_06_embedded_polarization_has_no_lead_diagonal(self):
        b=build_finite_system(Spec(length=25))
        self.assertEqual(np.linalg.norm(b['Pa'][4:,4:]),0.)
        self.assertEqual(np.linalg.norm(b['Pb'][4:,4:]),0.)

    def test_07_current_vertex_is_hermitian(self):
        b=build_finite_system(self.s)
        self.assertLess(np.linalg.norm(b['J0']-b['J0'].conj().T),1e-12)

    def test_08_electric_contact_is_hermitian(self):
        b=build_finite_system(self.s)
        self.assertLess(np.linalg.norm(b['JE']-b['JE'].conj().T),1e-12)

    def test_09_noncommuting_projected_positions(self):
        b=build_finite_system(self.s)
        self.assertGreater(np.linalg.norm(comm(b['Pa'][:4,:4],b['Pb'][:4,:4])),.1)
        self.assertGreater(np.linalg.norm(b['JE']),.01)

    def test_10_complex_finite_frequency_response(self):
        z=self.frames['symmetric']['sigma_device']
        self.assertTrue(np.isfinite([z.real,z.imag]).all())
        self.assertGreater(abs(z.imag),1e-4)

    def test_11_symmetric_landau_gauge_covariance(self):
        a=self.frames['symmetric']['sigma_device'];b=self.frames['landau_y']['sigma_device']
        self.assertLess(abs(a-b),2e-11)

    def test_12_noncommuting_orbital_frame_covariance(self):
        a=self.frames['symmetric']['sigma_device'];b=self.frames['orbital']['sigma_device']
        self.assertLess(abs(a-b),2e-11)

    def test_13_response_contains_electric_current_contact(self):
        r=self.frames['symmetric']
        self.assertGreater(abs(r['sigma_electric_contact']),1e-3)
        self.assertGreater(abs(r['sigma_device']-r['sigma_noncontact']),1e-3)

    def test_14_independent_lehmann_formula_is_finite(self):
        self.assertTrue(np.isfinite(self.frames['symmetric']['rhoE_eigenbasis']).all())

    def test_15_retarded_response_has_no_diagonal_population_drive(self):
        self.assertLess(np.linalg.norm(np.diag(self.frames['symmetric']['rhoE_eigenbasis'])),1e-14)

    def test_16_closed_system_total_ac_charge_zero(self):
        self.assertLess(abs(self.ward['total_charge_response']),1e-13)

    def test_17_device_charge_balanced_by_explicit_leads(self):
        self.assertGreater(abs(self.ward['device_charge_response']),.05)
        self.assertLess(abs(self.ward['device_charge_response']+self.ward['lead_charge_response']),1e-12)

    def test_18_exact_finite_lehmann_ward(self):
        self.assertLess(self.ward['ward_error'],1e-12)

    def test_19_charge_partition_identity(self):
        self.assertLess(self.ward['charge_partition_error'],1e-13)

    def test_20_ward_is_frame_covariant(self):
        for frame in ('landau_y','orbital'):
            w=lehmann_continuity_ward(Spec(length=40,eta=.03),frame)
            self.assertLess(w['ward_error'],1e-12)
            self.assertLess(abs(w['device_charge_response']-self.ward['device_charge_response']),2e-10)

    def test_21_static_isothermal_lehmann_matches_full_E_rediagonalization(self):
        self.assertLess(abs(self.static-self.static_fd),2e-9)

    def test_22_static_response_matches_frozen_reference(self):
        z=FROZEN['quick']['static_isothermal_sigma']
        self.assertLess(abs(self.static-complex(z['real'],z['imag'])),1e-9)

    def test_23_static_response_frame_covariance(self):
        for frame in ('landau_y','orbital'):
            x=static_isothermal_response(self.static_spec,frame)['sigma_device']
            self.assertLess(abs(x-self.static),2e-10)

    def test_24_static_dynamic_order_of_limits_is_not_silently_identified(self):
        dyn=lehmann_response(self.static_spec)['sigma_device']
        self.assertGreater(abs(dyn-self.static),1e-4)

    def test_25_static_current_FD_is_stable_to_step(self):
        a=static_finite_difference(self.static_spec,step=2e-4)
        b=static_finite_difference(self.static_spec,step=1e-4)
        self.assertLess(abs(a-b),2e-8)

    def test_26_static_mixed_EB_four_corner_symmetric(self):
        s=Spec(length=30,eta=.03,omega=0.)
        a=finite_difference_static_mixed(s,step=2e-3)
        b=finite_difference_static_mixed(s,'landau_y',step=2e-3)
        self.assertLess(abs(a-b),2e-8)

    def test_27_static_mixed_EB_orbital_frame(self):
        s=Spec(length=30,eta=.03,omega=0.)
        a=finite_difference_static_mixed(s,step=2e-3)
        b=finite_difference_static_mixed(s,'orbital',step=2e-3)
        self.assertLess(abs(a-b),2e-8)

    def test_28_full_B_five_point_derivative_is_stable(self):
        self.assertTrue(np.isfinite([self.long_B.real,self.long_B.imag]).all())

    def test_29_finite_lead_B_derivative_approaches_Keldysh(self):
        target=self.keldysh['mixed_dB_sigma_device']
        self.assertLess(abs(self.long_B-target),5e-4)

    def test_30_joint_length_broadening_sequence_improves_response(self):
        ref=self.keldysh['sigma_device']
        em=abs(self.medium['sigma_device']-ref)
        el=abs(self.long['sigma_device']-ref)
        ei=abs(self.intermediate['sigma_device']-ref)
        ep=abs(self.production['sigma_device']-ref)
        self.assertLess(el,em)
        self.assertLess(ei,el)
        self.assertLess(ep,ei)

    def test_31_medium_length_error_is_bounded(self):
        self.assertLess(abs(self.medium['sigma_device']-self.keldysh['sigma_device']),1.1e-3)

    def test_32_long_length_error_is_bounded(self):
        self.assertLess(abs(self.long['sigma_device']-self.keldysh['sigma_device']),9e-4)

    def test_33_production_lehmann_matches_continuum_within_resolution(self):
        self.assertLess(abs(self.production['sigma_device']-self.keldysh['sigma_device']),2.5e-4)

    def test_34_production_lehmann_regenerates_frozen_number(self):
        z=FROZEN['production']['sigma_device']
        val=self.production['sigma_device']
        self.assertLess(abs(val-complex(z['real'],z['imag'])),1e-7)

    def test_35_keldysh_reference_value_is_frozen(self):
        z=FROZEN['keldysh_reference']['sigma_device'];v=self.keldysh['sigma_device']
        self.assertLess(abs(v-complex(z['real'],z['imag'])),2e-9)

    def test_36_low_resolution_finite_lead_is_a_negative_control(self):
        low=lehmann_response(Spec(length=40,eta=.04))['sigma_device']
        self.assertGreater(abs(low-self.keldysh['sigma_device']),1e-3)

    def test_37_far_end_cap_effect_is_suppressed_by_long_leads(self):
        vals=[lehmann_response(Spec(length=120,eta=.02,end_cap=c))['sigma_device'] for c in (-1.1,0.,1.1)]
        self.assertLess(max(abs(v-vals[1]) for v in vals),1e-4)

    def test_38_omega_reflection_conjugation(self):
        pos=self.frames['symmetric']['sigma_device']
        neg=lehmann_response(replace(self.s,omega=-self.s.omega))['sigma_device']
        self.assertLess(abs(neg-pos.conjugate()),2e-10)

    def test_39_finite_response_changes_with_frequency(self):
        zero=lehmann_response(replace(self.s,omega=0.))['sigma_device']
        pos=self.frames['symmetric']['sigma_device']
        self.assertGreater(abs(zero-pos),1e-4)

    def test_40_invalid_lead_length_rejected(self):
        with self.assertRaises(ValueError): Spec(length=3)

    def test_41_invalid_eta_rejected(self):
        with self.assertRaises(ValueError): Spec(eta=0.)

    def test_42_invalid_temperature_rejected(self):
        with self.assertRaises(ValueError): Spec(temperature=0.)

    def test_43_invalid_static_frequency_rejected(self):
        with self.assertRaises(ValueError): static_isothermal_response(self.s)
        with self.assertRaises(ValueError): static_finite_difference(self.s)

    def test_44_invalid_finite_difference_rejected(self):
        with self.assertRaises(ValueError): finite_difference_B(self.s,step=0.)
        with self.assertRaises(ValueError): finite_difference_static_mixed(self.s,step=-1.)

    def test_45_status_keeps_original_paper_unresolved(self):
        self.assertIn('NOT_BULK',self.frames['symmetric']['status'])


if __name__=='__main__': unittest.main()
