"""v57 independent coherent 2/4-terminal NEGF and floating-voltage-probe checks."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.lead_transport_v57 import (
    LeadSpec,diamond_graph,_transmission_at_energy,_weights,transmissions,
    laplacian,four_terminal,device_response,full_finite_field_fd)
from qgmap.embedding_magnetic_v55 import peierls_field

class LandauerPhysicalV57(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pos,cls.h,cls.sites=diamond_graph(1)
        cls.spec=LeadSpec()
        cls.gammas=[.55]*4
        cls.a=device_response(cls.h,cls.pos,cls.sites,cls.gammas,cls.spec)
        cls.l=device_response(cls.h,cls.pos,cls.sites,cls.gammas,cls.spec,'landau_y')

    def test_01_true_4_terminal_transmissions_are_positive(self):
        T=self.a['T']
        self.assertGreater(np.min(T+np.eye(4)),0)
        self.assertLess(np.max(T),1.)
        np.testing.assert_allclose(np.diag(T),np.zeros(4),atol=0)

    def test_02_independent_caroli_matrix_trace_and_local_lead_formula(self):
        pb=peierls_field(self.h,self.pos,self.spec.field,self.spec.charge_over_hbar)
        E=.2;sites=self.sites;gamma=self.gammas
        T,dT=_transmission_at_energy(E,pb.hamiltonian,pb.magnetic_vertex,sites,gamma)
        Gamma=np.diag(gamma);gr=np.linalg.inv(E*np.eye(4)-pb.hamiltonian+.5j*Gamma)
        ga=gr.conj().T
        for alpha in range(4):
            for beta in range(4):
                if alpha==beta:continue
                P=np.diag([gamma[alpha] if i==sites[alpha] else 0 for i in range(4)])
                Q=np.diag([gamma[beta] if i==sites[beta] else 0 for i in range(4)])
                self.assertAlmostEqual(T[alpha,beta],np.trace(P@gr@Q@ga).real,delta=1e-13)
        self.assertGreater(np.max(abs(dT)),.05)

    def test_03_retarded_advanced_optical_theorem_current_conservation(self):
        p=peierls_field(self.h,self.pos,self.spec.field,self.spec.charge_over_hbar)
        E=.15;gam=np.diag(self.gammas)
        gr=np.linalg.inv(E*np.eye(4)-p.hamiltonian+.5j*gam)
        ga=gr.conj().T
        np.testing.assert_allclose(1j*(gr-ga),gr@gam@ga,atol=1e-13)
        np.testing.assert_allclose(self.a['T'].sum(axis=1),self.a['T'].sum(axis=0),atol=1e-12)

    def test_04_source_gauge_equivalence_of_every_transmission(self):
        np.testing.assert_allclose(self.a['T'],self.l['T'],atol=3e-13)
        np.testing.assert_allclose(self.a['dT_dB'],self.l['dT_dB'],atol=5e-12)

    def test_05_analytic_green_function_B_derivative_vs_full_five_point(self):
        fd=full_finite_field_fd(self.h,self.pos,self.sites,self.gammas,self.spec)
        np.testing.assert_allclose(fd['dT_dB'],self.a['dT_dB'],atol=1e-8)

    def test_06_independent_five_point_derivative_of_floating_hall_voltage(self):
        for gauge in ('symmetric','landau_y'):
            d=device_response(self.h,self.pos,self.sites,self.gammas,self.spec,gauge)
            fd=full_finite_field_fd(self.h,self.pos,self.sites,self.gammas,self.spec,gauge)
            self.assertAlmostEqual(fd['d_hall_ratio'],d['observables']['dHall_dB'],delta=1e-8)
            self.assertAlmostEqual(fd['d_longitudinal_conductance_q2_over_h'],d['observables']['dLongitudinal_dB'],delta=1e-8)

    def test_07_probe_currents_and_total_current_conservation(self):
        obs=self.a['observables']
        self.assertLess(obs['max_floating_probe_current'],1e-12)
        self.assertLess(obs['max_current_conservation_residual'],1e-12)
        self.assertLess(obs['max_derivative_current_conservation_residual'],1e-12)
        self.assertAlmostEqual(obs['lead_currents_h_over_q2'][0],-obs['lead_currents_h_over_q2'][2],delta=1e-12)

    def test_08_four_probe_magnetic_onsager_transmission_reciprocity(self):
        neg=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,field=-self.spec.field))
        np.testing.assert_allclose(self.a['T'],neg['T'].T,atol=2e-13)
        np.testing.assert_allclose(self.a['dT_dB'],-neg['dT_dB'].T,atol=2e-12)

    def test_09_two_terminal_magnetoconductance_even_field(self):
        sites=(0,2);gam=[.55,.55]
        p=device_response(self.h,self.pos,sites,gam,self.spec)
        n=device_response(self.h,self.pos,sites,gam,replace(self.spec,field=-self.spec.field))
        self.assertAlmostEqual(p['observables']['two_terminal_transmission_q2_over_h'],n['observables']['two_terminal_transmission_q2_over_h'],delta=2e-13)
        self.assertAlmostEqual(p['observables']['dTransmission_dB'],-n['observables']['dTransmission_dB'],delta=2e-12)

    def test_10_at_zero_magnetic_field_two_terminal_linear_slope_vanishes(self):
        z=device_response(self.h,self.pos,(0,2),[.55,.55],replace(self.spec,field=0))
        self.assertLess(abs(z['observables']['dTransmission_dB']),1e-12)

    def test_11_transverse_Hall_ratio_odd_for_mirror_symmetric_device(self):
        n=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,field=-self.spec.field))
        self.assertAlmostEqual(self.a['observables']['hall_ratio'],-n['observables']['hall_ratio'],delta=2e-13)
        self.assertAlmostEqual(self.a['observables']['dHall_dB'],n['observables']['dHall_dB'],delta=2e-12)
        self.assertGreater(abs(self.a['observables']['hall_ratio']),.05)

    def test_12_zero_field_4terminal_hall_voltage_vanishes(self):
        z=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,field=0))
        self.assertLess(abs(z['observables']['hall_ratio']),1e-12)
        self.assertGreater(abs(z['observables']['dHall_dB']),.01)

    def test_13_different_lead_coupling_changes_hall_and_does_not_need_gamma_bath(self):
        weak=device_response(self.h,self.pos,self.sites,[.2]*4,self.spec)
        strong=device_response(self.h,self.pos,self.sites,[1.1]*4,self.spec)
        self.assertGreater(abs(weak['observables']['dHall_dB']-strong['observables']['dHall_dB']),.04)
        self.assertNotIn('gamma',str(weak['limits']))

    def test_14_temp_smearing_changes_magnetic_hall_signal(self):
        cold=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,temperature=.04))
        warm=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,temperature=.35))
        self.assertGreater(abs(cold['observables']['hall_ratio']-warm['observables']['hall_ratio']),.06)

    def test_15_probes_are_open_and_reflection_consistent_under_origin_shift(self):
        for offset in ((2.,-.4),(-.5,1.6)):
            result=device_response(self.h,self.pos+offset,self.sites,self.gammas,self.spec)
            self.assertAlmostEqual(result['observables']['hall_ratio'],self.a['observables']['hall_ratio'],delta=1e-12)
            self.assertAlmostEqual(result['observables']['dHall_dB'],self.a['observables']['dHall_dB'],delta=2e-12)

    def test_16_thermodynamic_limit_not_faked_by_ring_size(self):
        signals=[]
        for n in (1,2,3,4):
            pos,h,sites=diamond_graph(n)
            a=device_response(h,pos,sites,[.55]*4,LeadSpec(energy_points=96))
            signals.append(a['observables']['dHall_dB'])
            self.assertLess(a['observables']['max_floating_probe_current'],1e-11)
        self.assertGreater(max(signals)-min(signals),5.)
        self.assertTrue(any(signals[i]*signals[i+1]<0 for i in range(3)))

    def test_17_energy_integration_normalized_and_convergent(self):
        energies,weights=_weights(self.spec)
        self.assertAlmostEqual(weights.sum(),np.tanh(self.spec.window_in_kT/2),delta=4e-11)
        a=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,energy_points=120))
        b=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,energy_points=220))
        self.assertAlmostEqual(a['observables']['hall_ratio'],b['observables']['hall_ratio'],delta=4e-6)

    def test_18_two_different_chemical_potentials_avoid_fake_quantization(self):
        a=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,mu=-.5))
        b=device_response(self.h,self.pos,self.sites,self.gammas,replace(self.spec,mu=.8))
        self.assertGreater(abs(a['observables']['hall_ratio']-b['observables']['hall_ratio']),.005)
        self.assertLess(abs(a['observables']['hall_ratio']),1)
        self.assertLess(abs(b['observables']['hall_ratio']),1)

    def test_19_reject_nonhermitian_invalid_leads_and_bad_temperature(self):
        with self.assertRaises(ValueError):diamond_graph(0)
        with self.assertRaises(ValueError):LeadSpec(temperature=0)
        with self.assertRaises(ValueError):LeadSpec(broadening=0)
        with self.assertRaises(ValueError):LeadSpec(energy_points=4)
        with self.assertRaises(ValueError):device_response(self.h,self.pos,self.sites,[.5,-.1,.5,.5],self.spec)
        with self.assertRaises(ValueError):device_response(self.h,self.pos,(0,0),[.5,.5],self.spec)
        with self.assertRaises(ValueError):device_response(self.h+1j*np.eye(4),self.pos,self.sites,self.gammas,self.spec)
        with self.assertRaises(ValueError):full_finite_field_fd(self.h,self.pos,self.sites,self.gammas,self.spec,step=0)
        with self.assertRaises(ValueError):four_terminal(np.eye(2),np.eye(2))

    def test_20_declared_device_not_bulk_intrinsic_or_v56_scalar_gamma(self):
        self.assertEqual(self.a['status'],'VERIFIED_COHERENT_LEAD_MODEL_NOT_BULK')
        self.assertEqual(self.a['source'],'wide_band_local_leads_finite_open_graph')
        self.assertNotIn('bulk',self.a['source'])
        self.assertGreater(self.a['observables']['longitudinal_conductance_q2_over_h'],0)

if __name__=='__main__':unittest.main()
