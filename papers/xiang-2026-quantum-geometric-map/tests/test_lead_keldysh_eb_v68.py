"""v68: noninteracting two-lead, energy-dependent finite-frequency Keldysh gates.

Tests verify a specific DEVICE polarization current operator in an explicitly
contacted model, not the bulk original-source E/B orbital Hall tensor.
"""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.lead_keldysh_eb_v68 import (
    Spec,comm,operator_data,semiinfinite_g,energy_grid,greens,response,
    keldysh_continuity,finite_difference_B,equilibrium_current,
    static_field_FD,static_mixed_FD)
from qgmap.dynamic_orbital_eb_v66 import dynamic_mixed_eb,benchmark_model

class V68NoninteractingKeldysh(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=Spec()
        cls.rs={f:response(cls.s,f) for f in ('symmetric','landau_y','orbital')}
        cls.w={f:keldysh_continuity(cls.s,f) for f in ('symmetric','landau_y','orbital')}
        cls.dc=replace(cls.s,omega=0.)
        cls.ref=cls.rs['symmetric']

    def test_01_exact_surface_quadratic(self):
        e=.13;t=1.1;g=semiinfinite_g(e,t,0)
        self.assertLess(abs(t*t*g*g-e*g+1),1e-14)

    def test_02_retarded_and_outside_band(self):
        self.assertLess(semiinfinite_g(.13,1.1,0).imag,0)
        self.assertAlmostEqual(semiinfinite_g(4.,1.1,0).imag,0)
        self.assertAlmostEqual(semiinfinite_g(-4.,1.1,0).imag,0)

    def test_03_full_projected_positions(self):
        o=self.ref['o'];self.assertGreater(np.linalg.norm(comm(o['Pa'],o['Pb'])),.35)

    def test_04_exact_green_optical_identity(self):
        o=operator_data(self.s)
        z=greens(.13,o,self.s)
        self.assertLess(np.linalg.norm(1j*(z['Gr']-z['Ga'])-z['A']),1e-9)

    def test_05_positive_spectral_weight(self):
        a=greens(.13,operator_data(self.s),self.s)['A']
        self.assertGreater(np.linalg.eigvalsh(a).min(),-1e-12)

    def test_06_total_spectral_normalization(self):
        o=operator_data(self.s);en,w=energy_grid(1.1,0,0,self.s.nodes_per_interval)
        avg=sum(weight*greens(e,o,self.s)['A'] for e,weight in zip(en,w))/(2*np.pi)
        self.assertLess(np.linalg.norm(avg-np.eye(4)),3e-8)

    def test_07_fermi_density_hermitian(self):
        self.assertLess(np.linalg.norm(self.ref['rho']-self.ref['rho'].conj().T),1e-12)

    def test_08_density_eigenvalues_physical(self):
        val=np.linalg.eigvalsh(self.ref['rho']);self.assertGreater(val.min(),0)
        self.assertLess(val.max(),1)

    def test_09_frame_density_covariance(self):
        r=self.ref
        for frame in ('landau_y','orbital'):
            x=self.rs[frame];U=x['o']['U'];self.assertLess(np.linalg.norm(x['rho']-U@r['rho']@U.conj().T),1e-11)

    def test_10_frame_ac_density_covariance(self):
        r=self.ref
        for frame in ('landau_y','orbital'):
            x=self.rs[frame];U=x['o']['U'];self.assertLess(np.linalg.norm(x['rhoE']-U@r['rhoE']@U.conj().T),1e-11)

    def test_11_frame_ac_response_covariance(self):
        for frame in ('landau_y','orbital'):
            self.assertLess(abs(self.rs[frame]['sigma_device']-self.ref['sigma_device']),1e-11)

    def test_12_frame_mixed_EB_response_covariance(self):
        for frame in ('landau_y','orbital'):
            self.assertLess(abs(self.rs[frame]['mixed_dB_sigma_device']-self.ref['mixed_dB_sigma_device']),1e-11)

    def test_13_ac_reference_converged(self):
        self.assertLess(abs(self.ref['sigma_device']-(-.0052926496159-.00119425739909j)),2e-9)

    def test_14_mixed_reference_converged(self):
        self.assertLess(abs(self.ref['mixed_dB_sigma_device']-(.00707005390013-.00414629084838j)),2e-9)

    def test_15_quadrature_order_convergence(self):
        co=response(replace(self.s,nodes_per_interval=512))
        self.assertLess(abs(co['sigma_device']-self.ref['sigma_device']),1e-9)
        self.assertLess(abs(co['mixed_dB_sigma_device']-self.ref['mixed_dB_sigma_device']),1e-9)

    def test_16_underresolved_resonance_negative_control(self):
        bad=response(replace(self.s,nodes_per_interval=64))
        self.assertGreater(abs(bad['sigma_device']-self.ref['sigma_device']),.001)

    def test_17_full_5point_B_derivative_symmetric(self):
        x=finite_difference_B(self.s,step=.005,points=5)
        self.assertLess(abs(x-self.ref['mixed_dB_sigma_device']),3e-10)

    def test_18_full_5point_B_derivative_landau(self):
        x=finite_difference_B(self.s,'landau_y',step=.005,points=5)
        self.assertLess(abs(x-self.ref['mixed_dB_sigma_device']),3e-10)

    def test_19_full_5point_B_derivative_orbital(self):
        x=finite_difference_B(self.s,'orbital',step=.005,points=5)
        self.assertLess(abs(x-self.ref['mixed_dB_sigma_device']),3e-10)

    def test_20_three_point_B_derivative_approaches_full(self):
        x=finite_difference_B(self.s,step=.002,points=3)
        self.assertLess(abs(x-self.ref['mixed_dB_sigma_device']),3e-8)

    def test_21_microscopic_two_energy_ward_all_frames(self):
        for f,v in self.w.items():self.assertLess(v['full_ward_error'],3e-9,f)

    def test_22_uncontacted_site_ac_continuity(self):
        for f,v in self.w.items():
            self.assertLess(abs(v['uncontacted_site_ward']),3e-9,f)
            self.assertLess(abs(v['uncontacted_site_collision']),1e-11,f)

    def test_23_finite_lead_collision_nonzero(self):
        self.assertGreater(self.w['symmetric']['lead_collision_norm'],.05)

    def test_24_global_charge_ward_with_lead_injection(self):
        v=self.w['symmetric']
        self.assertLess(abs(v['number_accumulation']-v['number_lead_injection']),3e-9)
        self.assertGreater(abs(v['number_lead_injection']),.1)

    def test_25_frequency_shifted_ward(self):
        v=keldysh_continuity(replace(self.s,omega=.13))
        self.assertLess(v['full_ward_error'],3e-9)

    def test_26_omega_reflection(self):
        neg=response(replace(self.s,omega=-.29))
        self.assertLess(abs(neg['sigma_device']-self.ref['sigma_device'].conjugate()),1e-8)

    def test_27_noncommuting_electric_current_contact(self):
        self.assertGreater(abs(self.ref['sigma_electric_contact']),.04)
        bad=response(self.s,include_electric_contact=False)
        self.assertGreater(abs(bad['sigma_device']-self.ref['sigma_device']),.04)

    def test_28_omit_magnetic_lead_contact_gauge_failure(self):
        bad=response(self.s,'orbital',omit_lead_magnetic_contacts=True)
        self.assertGreater(abs(bad['mixed_dB_sigma_device']-self.ref['mixed_dB_sigma_device']),.005)

    def test_29_omit_magnetic_position_contact_gauge_failure(self):
        bad=response(self.s,'orbital',include_all_magnetic_contacts=False)
        self.assertGreater(abs(bad['mixed_dB_sigma_device']-self.ref['mixed_dB_sigma_device']),.04)

    def test_30_full_magnetic_density_B_finite_difference(self):
        p=response(replace(self.s,field=self.s.field+.003))
        m=response(replace(self.s,field=self.s.field-.003))
        fd=(p['rho']-m['rho'])/.006
        self.assertLess(np.linalg.norm(fd-self.ref['rhoB']),4e-7)

    def test_31_full_magnetic_ac_density_B_finite_difference(self):
        p=response(replace(self.s,field=self.s.field+.003))
        m=response(replace(self.s,field=self.s.field-.003))
        fd=(p['rhoE']-m['rhoE'])/.006
        self.assertLess(np.linalg.norm(fd-self.ref['rhoEB']),3e-6)

    def test_32_independent_static_finite_E_full_current(self):
        d=response(self.dc)
        fd=static_field_FD(self.dc,step=.002)
        self.assertLess(abs(d['sigma_device']-fd),1e-8)

    def test_33_independent_static_finite_E_landau(self):
        d=response(self.dc,'landau_y')
        fd=static_field_FD(self.dc,'landau_y',step=.002)
        self.assertLess(abs(d['sigma_device']-fd),1e-8)

    def test_34_independent_static_finite_E_orbital(self):
        d=response(self.dc,'orbital')
        fd=static_field_FD(self.dc,'orbital',step=.002)
        self.assertLess(abs(d['sigma_device']-fd),1e-8)

    def test_35_independent_4corner_EB_DC_three_frames(self):
        for frame in ('symmetric','landau_y','orbital'):
            analytical=response(self.dc,frame)['mixed_dB_sigma_device']
            direct=static_mixed_FD(self.dc,frame,step=.002)
            self.assertLess(abs(analytical-direct),2e-8,frame)

    def test_36_explicit_lead_bath_not_scalar_gamma(self):
        h,p,px,py,k=benchmark_model()
        scalar=dynamic_mixed_eb(h,p,px,py)['sigma']
        self.assertGreater(abs(self.ref['sigma_device']-scalar),.005)
        self.assertIn('NONINTERACTING_LEAD_KELDYSH',self.ref['status'])

    def test_37_static_ward_has_lead_exchange_not_fictitious_zero(self):
        val=equilibrium_current(self.dc)
        self.assertGreater(abs(val),.006)
        w=keldysh_continuity(self.dc)
        self.assertLess(w['full_ward_error'],3e-9)

    def test_38_omega_changes_same_observable(self):
        at0=response(self.dc)['sigma_device']
        self.assertGreater(abs(at0-self.ref['sigma_device']),.008)

    def test_39_infinite_leads_not_WBA(self):
        o=operator_data(self.s)
        low=greens(-1.4,o,self.s)['Sigma']
        high=greens(.8,o,self.s)['Sigma']
        self.assertGreater(np.linalg.norm(low-high),.1)

    def test_40_fail_closed_invalid_inputs(self):
        for kw in ({'temperature':0},{'lead_hopping':0},{'nodes_per_interval':8},
                   {'omega':np.nan},{'mu':4},{'charge_over_hbar':0}):
            with self.assertRaises(ValueError):Spec(**kw)
        with self.assertRaises(ValueError):finite_difference_B(self.s,step=0)
        with self.assertRaises(ValueError):static_mixed_FD(self.s)

if __name__=='__main__':unittest.main()
