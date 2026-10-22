"""v67 physically different conserving semi-infinite lead scattering oracle."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.reservoir_ward_v67 import (
    ReservoirSpec,chain_surface,setup,reference,conductance,conductance_B,
    conductance_E,conductance_EB,local_polarization_bias,scattering,
    finite_difference,finite_difference_mixed)
from qgmap.dynamic_orbital_eb_v66 import benchmark_model
from qgmap.equilibrium_orbital_eb_v65 import operators,Spec

class ReservoirWardV67(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frames={f:reference(f) for f in ('symmetric','landau_y','orbital')}
        cls.base=cls.frames['symmetric']
        cls.sc={f:scattering(o) for f,o in cls.frames.items()}
        cls.local={f:local_polarization_bias(o) for f,o in cls.frames.items()}

    def test_01_exact_semi_infinite_surface_fixed_point(self):
        sp=ReservoirSpec(lead_regulator=.08)
        g=chain_surface(sp)
        z=sp.mu-sp.lead_onsite+1j*sp.lead_regulator
        self.assertLess(abs(g-1/(z-sp.lead_hopping**2*g)),1e-13)
        self.assertLess(g.imag,0)

    def test_02_semi_infinite_surface_independent_recursion(self):
        sp=ReservoirSpec(lead_regulator=.15)
        g=0j; z=sp.mu-sp.lead_onsite+1j*sp.lead_regulator
        for _ in range(600):g=1/(z-sp.lead_hopping**2*g)
        self.assertLess(abs(g-chain_surface(sp)),1e-11)

    def test_03_no_wide_band_constant_selfenergy(self):
        a=chain_surface(ReservoirSpec(mu=0.))
        b=chain_surface(ReservoirSpec(mu=.8))
        self.assertGreater(abs(a-b),.2)

    def test_04_positive_linewidth(self):
        for op in self.frames.values():
            for gamma in op['Gamma']:
                self.assertGreater(np.linalg.eigvalsh(gamma).min(),-1e-12)
                self.assertGreater(np.linalg.eigvalsh(gamma).max(),.01)

    def test_05_two_terminal_spectral_optical_identity(self):
        for op in self.frames.values():
            self.assertLess(op['spectral_identity_error'],1e-12)

    def test_06_S_unitarity(self):
        for s in self.sc.values():self.assertLess(s['unitarity_error'],1e-12)

    def test_07_independent_S_vs_Caroli(self):
        for name,op in self.frames.items():
            self.assertLess(abs(self.sc[name]['cross_probability']-conductance(op)),1e-12)

    def test_08_left_right_probability_ward(self):
        for sc in self.sc.values():
            self.assertLess(abs(sc['left_probability_sum']-1),1e-12)
            self.assertLess(abs(sc['right_probability_sum']-1),1e-12)

    def test_09_caroli_reciprocity(self):
        op=self.base;gl,gr=op['Gamma'];G=op['Gr'];A=op['Ga']
        self.assertLess(abs(np.trace(gl@G@gr@A)-np.trace(gr@G@gl@A)),1e-12)

    def test_10_two_terminal_T_frozen(self):
        self.assertLess(abs(conductance(self.base)-.9006260865577836),1e-12)

    def test_11_magnetic_gauge_T_invariant(self):
        for op in self.frames.values():
            self.assertLess(abs(conductance(op)-conductance(self.base)),1e-12)

    def test_12_magnetic_B_derivative_invariant(self):
        for op in self.frames.values():
            self.assertLess(abs(conductance_B(op)-conductance_B(self.base)),1e-12)

    def test_13_magnetic_EB_cross_invariant(self):
        for op in self.frames.values():
            self.assertLess(abs(conductance_EB(op)-conductance_EB(self.base)),1e-12)

    def test_14_B_derivative_independent_five_point(self):
        for frame,op in self.frames.items():
            self.assertLess(abs(conductance_B(op)-finite_difference(frame=frame)),2e-9)

    def test_15_E_derivative_independent_five_point(self):
        for frame,op in self.frames.items():
            self.assertLess(abs(conductance_E(op)-finite_difference(frame=frame,parameter='E')),2e-9)

    def test_16_EB_derivative_independent_sixteen_corner(self):
        for frame,op in self.frames.items():
            self.assertLess(abs(conductance_EB(op)-finite_difference_mixed(frame=frame)),3e-9)

    def test_17_symmetric_B_and_EB_numeric_reference(self):
        self.assertLess(abs(conductance_B(self.base)+.180487543398844),2e-12)
        self.assertLess(abs(conductance_EB(self.base)-.8840574786955058),2e-12)

    def test_18_complete_nonabelian_contact_magnetic_vertex(self):
        op=self.frames['orbital']
        self.assertGreater(np.linalg.norm(sum(op['Sd'])),.03)
        self.assertGreater(abs(conductance_B(op,omit_contact_B=True)-conductance_B(op)),.01)

    def test_19_complete_landau_magnetic_contact_vertex(self):
        op=self.frames['landau_y']
        self.assertGreater(np.linalg.norm(sum(op['Sd'])),.0001)
        self.assertGreater(abs(conductance_B(op,omit_contact_B=True)-conductance_B(op)),5e-5)

    def test_20_nonabelian_omit_lead_contact_fails_cross(self):
        op=self.frames['orbital']
        self.assertGreater(abs(conductance_EB(op,omit_magnetic_contacts=True)-conductance_EB(op)),.2)

    def test_21_landau_omit_lead_contact_fails_cross(self):
        op=self.frames['landau_y']
        self.assertGreater(abs(conductance_EB(op,omit_magnetic_contacts=True)-conductance_EB(op)),.005)

    def test_22_nonabelian_omit_position_contact_fails_cross(self):
        op=self.frames['orbital']
        self.assertGreater(abs(conductance_EB(op,omit_position_contact=True)-conductance_EB(op)),.09)

    def test_23_landau_omit_position_contact_fails_cross(self):
        op=self.frames['landau_y']
        self.assertGreater(abs(conductance_EB(op,omit_position_contact=True)-conductance_EB(op)),.015)

    def test_24_wrong_contact_transformation_fails_gauge_T(self):
        self.assertGreater(abs(conductance(reference('orbital',transform_contacts=False))-
                               conductance(self.base)),.006)
        self.assertGreater(abs(conductance(reference('landau_y',transform_contacts=False))-
                               conductance(self.base)),2e-5)

    def test_25_contact_selfenergy_frame_covariance(self):
        for frame in ('landau_y','orbital'):
            op=self.frames[frame];U=op['U']
            for i in range(2):
                self.assertLess(np.linalg.norm(op['S'][i]-U@self.base['S'][i]@U.conj().T),1e-12)

    def test_26_contact_selfenergy_B_derivative_direct_fd(self):
        for frame in ('landau_y','orbital'):
            op=self.frames[frame];h=1e-5
            spec=ReservoirSpec()
            p=reference(frame,spec=replace(spec,field=spec.field+h))
            m=reference(frame,spec=replace(spec,field=spec.field-h))
            for a in range(2):
                self.assertLess(np.linalg.norm(op['Sd'][a]-(p['S'][a]-m['S'][a])/(2*h)),1e-9)

    def test_27_uncontacted_site_charge_continuity(self):
        for item in self.local.values():self.assertLess(item['uncontacted_site3_continuity'],1e-12)

    def test_28_uncontacted_site_magnetic_charge_continuity(self):
        for item in self.local.values():self.assertLess(item['uncontacted_site3_continuity_B'],1e-12)

    def test_29_internal_Pa_current_is_distinct_from_terminal_current(self):
        self.assertGreater(abs(self.local['symmetric']['internal_pa_bias_current']-conductance(self.base)),1)

    def test_30_internal_Pa_bias_current_frame_covariance(self):
        ref=self.local['symmetric']
        for item in self.local.values():
            self.assertLess(abs(item['internal_pa_bias_current']-ref['internal_pa_bias_current']),1e-12)
            self.assertLess(abs(item['internal_pa_bias_current_B']-ref['internal_pa_bias_current_B']),1e-12)

    def test_31_internal_Pa_current_B_finite_diff(self):
        s=ReservoirSpec();h=1e-5
        for frame,ref in self.local.items():
            p=local_polarization_bias(reference(frame,spec=replace(s,field=s.field+h)))
            m=local_polarization_bias(reference(frame,spec=replace(s,field=s.field-h)))
            val=(p['internal_pa_bias_current']-m['internal_pa_bias_current'])/(2*h)
            self.assertLess(abs(val-ref['internal_pa_bias_current_B']),2e-9)

    def test_32_projection_physical_Pa_Pb_same_inputs_as_v66(self):
        h,xy,pa,pb,k=benchmark_model()
        sp=Spec(mu=.13,temperature=.37,field=.23,charge_over_hbar=-.84)
        opa=operators(h,xy,pa,sp)
        opb=operators(h,xy,pb,sp)
        self.assertLess(np.linalg.norm(self.base['Pa']-opa['P']),1e-12)
        self.assertLess(np.linalg.norm(self.base['Pb']-opb['P']),1e-12)

    def test_33_gate_origin_shift_requires_reservoir_reference_shift(self):
        h,xy,pa,pb,k=benchmark_model();s=ReservoirSpec();c=1.23
        shifted=replace(s,mu=s.mu-s.electric_gate*c,lead_onsite=s.lead_onsite-s.electric_gate*c)
        val=setup(h,xy,pa,pb+c*np.eye(4),shifted)
        self.assertLess(abs(conductance(val)-conductance(self.base)),1e-12)
        wrong=setup(h,xy,pa,pb+c*np.eye(4),s)
        self.assertGreater(abs(conductance(wrong)-conductance(self.base)),.01)

    def test_34_intermediate_energy_sweep_analytic_surface(self):
        vals=[conductance(reference(spec=replace(ReservoirSpec(),mu=m))) for m in (-.7,0.,.13,.7)]
        self.assertTrue(np.isfinite(vals).all())
        self.assertGreater(max(vals)-min(vals),.03)

    def test_35_regulator_control_not_gamma_bath(self):
        s=ReservoirSpec(lead_regulator=0)
        op=reference(spec=s)
        self.assertLess(scattering(op)['unitarity_error'],1e-12)
        self.assertGreater(conductance(op),.7)

    def test_36_no_channel_outside_lead_band(self):
        with self.assertRaises(ValueError):ReservoirSpec(mu=2.3)

    def test_37_fail_negative_hopping(self):
        with self.assertRaises(ValueError):ReservoirSpec(lead_hopping=0)

    def test_38_fail_negative_lead_regulator(self):
        with self.assertRaises(ValueError):ReservoirSpec(lead_regulator=-1e-3)

    def test_39_fail_nonfinite_parameters(self):
        with self.assertRaises(ValueError):ReservoirSpec(field=float('nan'))

    def test_40_fail_fd_invalid_args(self):
        with self.assertRaises(ValueError):finite_difference(points=4)
        with self.assertRaises(ValueError):finite_difference(parameter='mu')
        with self.assertRaises(ValueError):finite_difference_mixed(step=0)
        with self.assertRaises(ValueError):finite_difference_mixed(points=4)

    def test_41_proof_scope_not_eb_equivalence(self):
        # The observable is a two-terminal reservoir B/gate derivative, not sigma_xy.
        self.assertNotAlmostEqual(conductance_EB(self.base),-.006402948002937392)

if __name__=='__main__':unittest.main()
