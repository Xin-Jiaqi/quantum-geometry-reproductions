"""v71 source-defined spin-index and finite-lead observable separation tests."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.source_observable_bridge_v71 import (
    Spec,finite_contact_identity,equilibrium_full_current,source_eq20_eq31_sign,
    pauli_helicity_gate,capability_gate,REQUIRED_BY_CHANNEL)
from qgmap.finite_lead_lehmann_v69 import build_finite_system,lehmann_response

class V71SourceObservableBridge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec(length=22,eta=.025,omega=.29)
        cls.frames={frame:finite_contact_identity(cls.s,frame) for frame in ('symmetric','landau_y','orbital')}
        rng=np.random.default_rng(734)
        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        b=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        cls.r=(a+a.conj().T)/2
        cls.spin=(b+b.conj().T)/2
        cls.e=np.array([-.8,-.21,.51,1.2]);cls.df=np.array([.3,-.11,.19,-.25])
        cls.sign=source_eq20_eq31_sign(cls.e,cls.r,cls.spin,cls.df)

    def test_01_finite_contact_nonzero(self):
        self.assertGreater(self.frames['symmetric']['contact_operator_norm'],.7)

    def test_02_exact_current_decomposition(self):
        for r in self.frames.values():self.assertLess(r['current_decomposition_residual'],1e-13)

    def test_03_dynamic_response_decomposition(self):
        for r in self.frames.values():
            self.assertLess(abs(r['dynamic_full']-r['dynamic_device']-r['dynamic_contact']),2e-14)

    def test_04_full_dynamic_polarization_ward(self):
        for r in self.frames.values():self.assertLess(r['dynamic_ward_residual'],2e-13)

    def test_05_static_isothermal_complete_current_zero(self):
        for r in self.frames.values():self.assertLess(abs(r['static_isothermal_full']),2e-12)

    def test_06_static_device_current_nonzero(self):
        self.assertGreater(abs(self.frames['symmetric']['static_isothermal_device']),.002)

    def test_07_static_contact_exactly_compensates(self):
        r=self.frames['symmetric']
        self.assertLess(abs(r['static_isothermal_device']+r['static_isothermal_contact']),2e-12)

    def test_08_static_reference_device(self):
        self.assertAlmostEqual(self.frames['symmetric']['static_isothermal_device'].real,.00299637669423265,places=11)

    def test_09_static_contact_reference(self):
        self.assertAlmostEqual(self.frames['symmetric']['static_isothermal_contact'].real,-.00299637669423269,places=11)

    def test_10_full_at_finite_E_zero(self):
        for E in (-.07,-.003,0,.002,.09):
            r=equilibrium_full_current(replace(self.s,omega=0),E=E)
            self.assertLess(abs(r['full']),3e-12)

    def test_11_full_at_finite_B_zero(self):
        for b in (-.31,0.,.14,.53):
            r=equilibrium_full_current(replace(self.s,omega=0,field=b),E=.025)
            self.assertLess(abs(r['full']),3e-12)

    def test_12_full_static_finite_E_derivative_zero(self):
        spec=replace(self.s,omega=0)
        h=3e-4
        p=equilibrium_full_current(spec,E=h)['full'];m=equilibrium_full_current(spec,E=-h)['full']
        self.assertLess(abs((p-m)/(2*h)),3e-10)

    def test_13_device_static_finite_E_matches_analytic(self):
        spec=replace(self.s,omega=0);h=3e-4
        p=equilibrium_full_current(spec,E=h)['device'];m=equilibrium_full_current(spec,E=-h)['device']
        actual=(p-m)/(2*h)
        self.assertLess(abs(actual-self.frames['symmetric']['static_isothermal_device']),5e-9)

    def test_14_contact_static_finite_E_matches_analytic(self):
        spec=replace(self.s,omega=0);h=3e-4
        p=equilibrium_full_current(spec,E=h)['contact'];m=equilibrium_full_current(spec,E=-h)['contact']
        actual=(p-m)/(2*h)
        self.assertLess(abs(actual-self.frames['symmetric']['static_isothermal_contact']),5e-9)

    def test_15_equilibrium_current_partition_at_finite_E(self):
        for E in (-.03,.02):
            r=equilibrium_full_current(replace(self.s,omega=0),E=E)
            self.assertLess(abs(r['device']+r['contact']-r['full']),2e-13)

    def test_16_static_B_derivative_of_complete_current_zero(self):
        h=1e-3;spec=replace(self.s,omega=0)
        def val(B):return finite_contact_identity(replace(spec,field=B))['static_isothermal_full']
        self.assertLess(abs((val(.23+h)-val(.23-h))/(2*h)),1e-9)

    def test_17_dynamic_all_frames_equal(self):
        a=self.frames['symmetric']['dynamic_full']
        for r in self.frames.values():self.assertLess(abs(a-r['dynamic_full']),2e-12)

    def test_18_dynamic_contact_all_frames_equal(self):
        a=self.frames['symmetric']['dynamic_contact']
        for r in self.frames.values():self.assertLess(abs(a-r['dynamic_contact']),2e-12)

    def test_19_dynamic_device_reproduces_v69(self):
        ref=lehmann_response(self.s)['sigma_device']
        self.assertLess(abs(self.frames['symmetric']['dynamic_device']-ref),2e-13)

    def test_20_dynamic_full_differs_from_device(self):
        r=self.frames['symmetric'];self.assertGreater(abs(r['dynamic_full']-r['dynamic_device']),.02)

    def test_21_electric_current_contact_essential(self):
        r=self.frames['symmetric'];self.assertGreater(abs(r['electric_current_contact']),.05)

    def test_22_dynamic_zero_frequency_eta_resolution_nonzero(self):
        r=finite_contact_identity(replace(self.s,omega=0))
        self.assertGreater(abs(r['dynamic_full']),1e-5)
        self.assertLess(r['dynamic_ward_residual'],1e-12)

    def test_23_zero_frequency_eta_suppresses_full(self):
        a=finite_contact_identity(replace(self.s,omega=0,eta=.02))
        b=finite_contact_identity(replace(self.s,omega=0,eta=.01))
        self.assertLess(abs(b['dynamic_full']),abs(a['dynamic_full']))

    def test_24_literal_eq20_inferred_Z_identity(self):
        self.assertLess(self.sign['algebra_residual'],2e-13)

    def test_25_printed_eq31_disagrees_in_generic_spin_model(self):
        self.assertGreater(abs(self.sign['sign_discrepancy']),.2)

    def test_26_printed_and_inferred_branches_opposite(self):
        self.assertLess(abs(self.sign['printed_eq31']+self.sign['inferred_Z']),1e-14)

    def test_27_random_U1_gauge_covariance(self):
        phases=np.exp(1j*np.array([.12,-.41,.86,-.21]))
        U=np.diag(phases)
        r=U.conj().T@self.r@U;s=U.conj().T@self.spin@U
        a=source_eq20_eq31_sign(self.e,r,s,self.df)
        self.assertLess(abs(a['direct_eq20']-self.sign['direct_eq20']),1e-13)
        self.assertLess(abs(a['printed_eq31']-self.sign['printed_eq31']),1e-13)

    def test_28_permuted_band_labels_invariant(self):
        perm=np.array([2,0,3,1])
        a=source_eq20_eq31_sign(self.e[perm],self.r[np.ix_(perm,perm)],self.spin[np.ix_(perm,perm)],self.df[perm])
        self.assertLess(abs(a['direct_eq20']-self.sign['direct_eq20']),1e-13)

    def test_29_constant_derivative_control_gives_zero(self):
        a=source_eq20_eq31_sign(self.e,self.r,self.spin,np.full(4,.3))
        self.assertLess(abs(a['direct_eq20']),1e-13)

    def test_30_spin_zero_control(self):
        a=source_eq20_eq31_sign(self.e,self.r,np.zeros_like(self.spin),self.df)
        self.assertEqual(a['direct_eq20'],0j)

    def test_31_offdiagonal_position_zero_control(self):
        a=source_eq20_eq31_sign(self.e,np.diag(np.diag(self.r)),self.spin,self.df)
        self.assertEqual(a['direct_eq20'],0j)

    def test_32_original_main_sm_H_unitary(self):
        self.assertLess(pauli_helicity_gate()['hamiltonian_equivalence_error'],1e-13)

    def test_33_spin_must_transform_under_helicity_change(self):
        p=pauli_helicity_gate();self.assertLess(p['transformed_spin_error'],1e-13)
        self.assertGreater(p['holding_spin_fixed_error'],2)

    def test_34_sign_fail_closed_degeneracy(self):
        with self.assertRaises(ValueError):source_eq20_eq31_sign([1,1],np.eye(2),np.eye(2),[.1,.2])

    def test_35_sign_fail_closed_nonhermitian(self):
        with self.assertRaises(ValueError):source_eq20_eq31_sign(self.e,self.r+1j*np.eye(4),self.spin,self.df)

    def test_36_no_spin_no_spin_PHE(self):
        with self.assertRaisesRegex(ValueError,'spin_zeeman'):
            capability_gate('spin_phe',{'H_k':'given','velocity_k':'given'})

    def test_37_no_full_wannier_no_orbital_MNHE(self):
        with self.assertRaisesRegex(ValueError,'full_wannier_position'):
            capability_gate('orbital_mnhe',{'H_k':'given','velocity_k':'given'})

    def test_38_orbital_embedding_explicit_verification_required(self):
        supplied={k:'given' for k in REQUIRED_BY_CHANNEL['orbital_mnhe']}
        with self.assertRaisesRegex(ValueError,'verified_full_magnetic_Wannier_embedding'):
            capability_gate('orbital_mnhe',supplied)

    def test_39_even_complete_inventory_cannot_certify_bulk(self):
        supplied={k:'given' for k in REQUIRED_BY_CHANNEL['orbital_mnhe']}
        r=capability_gate('orbital_mnhe',supplied,complete_embedding=True)
        self.assertFalse(r['ready_for_author_claim']);self.assertFalse(r['proof_of_bulk_mapping'])

    def test_40_unknown_channel_rejected(self):
        with self.assertRaises(ValueError):capability_gate('all_bulk',{})

    def test_41_none_field_rejected(self):
        p={k:'given' for k in REQUIRED_BY_CHANNEL['spin_phe']};p['spin_zeeman']=None
        with self.assertRaises(ValueError):capability_gate('spin_phe',p)

    def test_42_module_status_not_original_source(self):
        self.assertIn('NOT_BULK_TRANSPORT',self.frames['symmetric']['status'])
        self.assertIn('NOT_AUTHOR_ERRATUM',self.sign['status'])

if __name__=='__main__':unittest.main()
