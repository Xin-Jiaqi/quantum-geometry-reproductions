"""v66 independent finite AC EB projected position gate; physical bath still CONDITIONAL."""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.dynamic_orbital_eb_v66 import (Spec,benchmark_model,ac_state,dynamic_mixed_eb,
    full_field_finite_difference,ac_spectral_independent,independent_dc_steady_current,
    independent_dc_mixed_fd,comm)
from qgmap.equilibrium_orbital_eb_v65 import (fermi,finite_fermi_frechet)


class V66ProjectedDynamic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.H,cls.xy,cls.Pa,cls.Pb,cls.K=benchmark_model()
        cls.s=Spec()
        cls.r=dynamic_mixed_eb(cls.H,cls.xy,cls.Pa,cls.Pb)
        cls.landau=dynamic_mixed_eb(cls.H,cls.xy,cls.Pa,cls.Pb,frame='landau_y')
        cls.orbital=dynamic_mixed_eb(cls.H,cls.xy,cls.Pa,cls.Pb,frame='orbital',K=cls.K)

    def test_01_projected_positions_noncommute(self):
        self.assertGreater(np.linalg.norm(comm(self.Pa,self.Pb)),.1)

    def test_02_defined_complex_response(self):
        self.assertLess(abs(self.r['sigma']-(-.0002810991670626692-.015989993886819445j)),1e-12)

    def test_03_signed_dynamic_mixed(self):
        self.assertLess(abs(self.r['mixed_dsigma_dB']-(-.006402948002937392+.004503188579708714j)),1e-12)

    def test_04_landau_dynamic_frame(self):
        self.assertLess(abs(self.r['mixed_dsigma_dB']-self.landau['mixed_dsigma_dB']),1e-12)
        self.assertLess(abs(self.r['sigma']-self.landau['sigma']),1e-12)

    def test_05_noncommuting_orbital_frame(self):
        self.assertLess(abs(self.r['mixed_dsigma_dB']-self.orbital['mixed_dsigma_dB']),1e-12)
        self.assertLess(abs(self.r['sigma']-self.orbital['sigma']),1e-12)

    def test_06_independent_eigenbasis_double_sum(self):
        for frame,k in [('symmetric',None),('landau_y',None),('orbital',self.K)]:
            x=ac_spectral_independent(self.H,self.xy,self.Pa,self.Pb,self.s,frame,k)
            ref=ac_state(self.H,self.xy,self.Pa,self.Pb,self.s,frame,k)
            self.assertLess(abs(x-ref['sigma']),1e-12)

    def test_07_independent_5pt_magnetic_in_symmetric(self):
        x=full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,step=1e-4)
        self.assertLess(abs(x-self.r['mixed_dsigma_dB']),2e-10)

    def test_08_independent_5pt_landau(self):
        x=full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,frame='landau_y')
        self.assertLess(abs(x-self.r['mixed_dsigma_dB']),2e-10)

    def test_09_independent_5pt_orbital(self):
        x=full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,frame='orbital',K=self.K)
        self.assertLess(abs(x-self.r['mixed_dsigma_dB']),2e-10)

    def test_10_three_vs_five_point(self):
        x=full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,step=1e-4,points=3)
        self.assertLess(abs(x-self.r['mixed_dsigma_dB']),2e-8)

    def test_11_seven_terms_sum(self):
        self.assertEqual(len(self.r['components']),7)
        self.assertLess(abs(sum(self.r['components'].values())-self.r['mixed_dsigma_dB']),1e-14)

    def test_12_source_position_contact_finite(self):
        self.assertGreater(abs(self.landau['components']['source_position_dPb']),.002)

    def test_13_current_position_contact_finite(self):
        self.assertGreater(abs(self.landau['components']['current_position_dPa']),.02)

    def test_14_electric_current_contact_B_finite(self):
        self.assertGreater(abs(self.landau['components']['electric_current_contact_dB']),.02)

    def test_15_electric_current_contact_density_finite(self):
        self.assertGreater(abs(self.r['components']['electric_current_contact_rhoB']),.04)

    def test_16_omit_position_source_fails_FD(self):
        self.assertGreater(abs(self.landau['mixed_dsigma_dB']-self.landau['components']['source_position_dPb']-
                               full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,frame='landau_y')),.001)

    def test_17_omit_current_position_fails_FD(self):
        self.assertGreater(abs(self.landau['mixed_dsigma_dB']-self.landau['components']['current_position_dPa']-
                               full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,frame='landau_y')),.02)

    def test_18_omit_electric_contact_derivative_fails(self):
        self.assertGreater(abs(self.landau['mixed_dsigma_dB']-self.landau['components']['electric_current_contact_dB']-
                               full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,frame='landau_y')),.02)

    def test_19_omit_current_E_contact_fails_physics(self):
        ref=ac_state(self.H,self.xy,self.Pa,self.Pb)
        wrong=ac_state(self.H,self.xy,self.Pa,self.Pb,electric_contact=False)
        self.assertGreater(abs(ref['sigma']-wrong['sigma']),.05)

    def test_20_rho_E_liouville_residual(self):
        s=ac_state(self.H,self.xy,self.Pa,self.Pb)
        residual=s['z']*s['rhoE']+1j*comm(s['H'],s['rhoE'])-s['source']
        self.assertLess(np.linalg.norm(residual),1e-12)

    def test_21_complex_frequency_reality(self):
        opposite=ac_state(self.H,self.xy,self.Pa,self.Pb,replace(self.s,omega=-self.s.omega))['sigma']
        self.assertLess(abs(opposite-self.r['sigma'].conjugate()),1e-12)

    def test_22_mixed_frequency_reality(self):
        opposite=dynamic_mixed_eb(self.H,self.xy,self.Pa,self.Pb,replace(self.s,omega=-self.s.omega))['mixed_dsigma_dB']
        self.assertLess(abs(opposite-self.r['mixed_dsigma_dB'].conjugate()),1e-12)

    def test_23_independent_superoperator_DC_sigma(self):
        spec=replace(self.s,omega=0)
        h=2e-5
        i1=independent_dc_steady_current(self.H,self.xy,self.Pa,self.Pb,spec,E=h)
        i2=independent_dc_steady_current(self.H,self.xy,self.Pa,self.Pb,spec,E=-h)
        sigma=(i1-i2)/(2*h)
        self.assertLess(abs(sigma-ac_state(self.H,self.xy,self.Pa,self.Pb,spec)['sigma']),1e-9)

    def test_24_DC_four_corner_reference(self):
        spec=replace(self.s,omega=0)
        a=dynamic_mixed_eb(self.H,self.xy,self.Pa,self.Pb,spec)['mixed_dsigma_dB']
        b=independent_dc_mixed_fd(self.H,self.xy,self.Pa,self.Pb,spec)
        self.assertLess(abs(a-b),2e-8)

    def test_25_DC_landau_four_corner_reference(self):
        spec=replace(self.s,omega=0)
        a=dynamic_mixed_eb(self.H,self.xy,self.Pa,self.Pb,spec,'landau_y')['mixed_dsigma_dB']
        b=independent_dc_mixed_fd(self.H,self.xy,self.Pa,self.Pb,spec,'landau_y')
        self.assertLess(abs(a-b),2e-8)

    def test_26_DC_orbital_four_corner_reference(self):
        spec=replace(self.s,omega=0)
        a=dynamic_mixed_eb(self.H,self.xy,self.Pa,self.Pb,spec,'orbital',self.K)['mixed_dsigma_dB']
        b=independent_dc_mixed_fd(self.H,self.xy,self.Pa,self.Pb,spec,'orbital',self.K)
        self.assertLess(abs(a-b),2e-8)

    def test_27_static_gamma_to_zero_current_cancelled(self):
        s=replace(self.s,omega=0,gamma=1e-6)
        v=ac_state(self.H,self.xy,self.Pa,self.Pb,s)['sigma']
        wrong=ac_state(self.H,self.xy,self.Pa,self.Pb,s,electric_contact=False)['sigma']
        self.assertLess(abs(v),5e-7)
        self.assertGreater(abs(wrong),.05)

    def test_28_static_gamma_to_zero_mixed_cancelled(self):
        s=replace(self.s,omega=0,gamma=1e-5)
        self.assertLess(abs(dynamic_mixed_eb(self.H,self.xy,self.Pa,self.Pb,s)['mixed_dsigma_dB']),1e-6)

    def test_29_spinless_TR_zero_B_onsager(self):
        H,xy,Pa,Pb,_=benchmark_model(False)
        s=replace(self.s,field=0)
        x=ac_state(H,xy,Pa,Pb,s)['sigma']
        y=ac_state(H,xy,Pb,Pa,s)['sigma']
        self.assertLess(abs(x-y),1e-12)

    def test_30_spinless_TR_onsager_B_reversal(self):
        H,xy,Pa,Pb,_=benchmark_model(False)
        x=ac_state(H,xy,Pa,Pb,replace(self.s,field=.23))['sigma']
        y=ac_state(H,xy,Pb,Pa,replace(self.s,field=-.23))['sigma']
        self.assertLess(abs(x-y),1e-12)

    def test_31_real_position_omits_E_current_contact(self):
        H,xy,_,_,_=benchmark_model(False)
        pa=np.diag(self.s.charge_over_hbar*xy[:,0]);pb=np.diag(self.s.charge_over_hbar*xy[:,1])
        x=ac_state(H,xy,pa,pb)
        self.assertLess(abs(np.trace(x['JE']@x['rho'])),1e-13)

    def test_32_under_specified_gamma_not_infer_conserving(self):
        self.assertIn('CONDITIONAL',self.r['status'])
        self.assertIn('NOT_CONSERVING_STREDA',self.r['status'])

    def test_33_invalid_gamma_refusal(self):
        with self.assertRaises(ValueError):Spec(gamma=0)

    def test_34_invalid_temperature_refusal(self):
        with self.assertRaises(ValueError):Spec(temperature=0)

    def test_35_invalid_mag_frame_refusal(self):
        with self.assertRaises(ValueError):ac_state(self.H,self.xy,self.Pa,self.Pb,frame='not_a_frame')

    def test_36_invalid_five_point_refusal(self):
        with self.assertRaises(ValueError):full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,step=-1)
        with self.assertRaises(ValueError):full_field_finite_difference(self.H,self.xy,self.Pa,self.Pb,points=4)

    def test_37_independent_dc_mixed_refuses_ac(self):
        with self.assertRaises(ValueError):independent_dc_mixed_fd(self.H,self.xy,self.Pa,self.Pb,self.s)

    def test_38_DC_superoperator_hermiticity(self):
        spec=replace(self.s,omega=0)
        # Real steady current from a static electric field, no fictitious imaginary current.
        val=independent_dc_steady_current(self.H,self.xy,self.Pa,self.Pb,spec,E=.02)
        self.assertLess(abs(val.imag),1e-12)

    def test_39_offdiag_fermi_B_nondegenerate(self):
        s=ac_state(self.H,self.xy,self.Pa,self.Pb)
        M=s['H']*0+np.eye(4)
        f=finite_fermi_frechet(s['H'],M,self.s.mu,self.s.temperature)
        eig,U=np.linalg.eigh(s['H']);p=fermi(eig,self.s.mu,self.s.temperature)
        exact=(U*np.array(-p*(1-p)/self.s.temperature)[None,:])@U.conj().T
        self.assertLess(np.linalg.norm(f-exact),1e-12)

    def test_40_negative_control_contact_magnetic_dPa(self):
        # Two basis frames redistribute each of the 7 diagram groups.
        self.assertGreater(abs(self.r['components']['current_magnetic_H']-
                               self.orbital['components']['current_magnetic_H']),.15)


if __name__=='__main__': unittest.main()
