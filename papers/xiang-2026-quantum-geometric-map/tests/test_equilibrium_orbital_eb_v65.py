"""v65 source-independent equilibrium Maxwell, projected position and ensemble gates.

40 named tests. No assertion of real-material/original-paper EB conductivity.
"""
import unittest
from dataclasses import replace
import numpy as np
from qgmap.equilibrium_orbital_eb_v65 import (
    Spec,benchmark_graph,full_projected_polarization,operators,fermi,
    finite_fermi_frechet,grand_potential,equilibrium_derivatives,mixed_static_EB,
    finite_difference_mixed,canonical_mixed_static_EB,finite_difference_canonical)
from qgmap.embedding_magnetic_v55 import peierls_field


class V65EquilibriumOrbitalEB(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.H,cls.xy,cls.D=benchmark_graph()
        cls.s=Spec()
        cls.P=full_projected_polarization(cls.xy,cls.s.charge_over_hbar,'y',cls.D)
        cls.K=np.array([[.3,.21+.1j,0,0],[.21-.1j,-.4,.11,0],
                        [0,.11,.2,.17+.08j],[0,0,.17-.08j,-.1]],complex)
        cls.ref=mixed_static_EB(cls.H,cls.xy,cls.P)
        cls.canon=canonical_mixed_static_EB(cls.H,cls.xy,cls.P)

    def test_01_peierls_exact(self):
        a=operators(self.H,self.xy,self.P)
        b=peierls_field(self.H,self.xy,self.s.field,self.s.charge_over_hbar)
        self.assertLess(np.linalg.norm(a['H']-b.hamiltonian),1e-13)
        self.assertLess(np.linalg.norm(a['M']-b.magnetic_vertex),1e-13)

    def test_02_offdiagonal_position(self):
        self.assertAlmostEqual(abs(self.P[0,1]),abs(self.D[0,1]),places=12)
        self.assertGreater(abs(self.P[2,3]),.15)

    def test_03_offdiagonal_response_not_center_only(self):
        bare=full_projected_polarization(self.xy,self.s.charge_over_hbar)
        m=mixed_static_EB(self.H,self.xy,bare)
        self.assertGreater(abs(m['mixed_d2F_dEdB']-self.ref['mixed_d2F_dEdB']),.001)

    def test_04_landau_operator_transformation(self):
        sy=operators(self.H,self.xy,self.P)
        ld=operators(self.H,self.xy,self.P,frame='landau_y');U=ld['U']
        self.assertLess(np.linalg.norm(ld['H']-U@sy['H']@U.conj().T),1e-13)
        self.assertLess(np.linalg.norm(ld['P']-U@sy['P']@U.conj().T),1e-13)

    def test_05_landau_dH_full_B_difference(self):
        d=1e-5;a=operators(self.H,self.xy,self.P,frame='landau_y')
        p=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field+d),frame='landau_y')
        m=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field-d),frame='landau_y')
        self.assertLess(np.linalg.norm(a['M']-(p['H']-m['H'])/(2*d)),1e-9)

    def test_06_landau_dP_full_B_difference(self):
        d=1e-5;a=operators(self.H,self.xy,self.P,frame='landau_y')
        p=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field+d),frame='landau_y')
        m=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field-d),frame='landau_y')
        self.assertLess(np.linalg.norm(a['dP']-(p['P']-m['P'])/(2*d)),1e-9)

    def test_07_landau_full_gauge_maxwell(self):
        ld=mixed_static_EB(self.H,self.xy,self.P,frame='landau_y')
        self.assertLess(abs(ld['mixed_d2F_dEdB']-self.ref['mixed_d2F_dEdB']),1e-12)

    def test_08_landau_omit_dP_fails(self):
        ld=mixed_static_EB(self.H,self.xy,self.P,frame='landau_y')
        self.assertGreater(abs(ld['fermi_density_derivative']-self.ref['mixed_d2F_dEdB']),.025)
        self.assertGreater(abs(ld['position_contact_dP_dB']),.025)

    def test_09_noncommuting_B_frame_unitary(self):
        a=operators(self.H,self.xy,self.P,frame='orbital',K=self.K)
        self.assertLess(np.linalg.norm(a['U'].conj().T@a['U']-np.eye(4)),1e-12)

    def test_10_nonabelian_H_derivative(self):
        d=1e-5;a=operators(self.H,self.xy,self.P,frame='orbital',K=self.K)
        p=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field+d),frame='orbital',K=self.K)
        m=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field-d),frame='orbital',K=self.K)
        self.assertLess(np.linalg.norm(a['M']-(p['H']-m['H'])/(2*d)),1e-9)

    def test_11_nonabelian_P_derivative(self):
        d=1e-5;a=operators(self.H,self.xy,self.P,frame='orbital',K=self.K)
        p=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field+d),frame='orbital',K=self.K)
        m=operators(self.H,self.xy,self.P,replace(self.s,field=self.s.field-d),frame='orbital',K=self.K)
        self.assertLess(np.linalg.norm(a['dP']-(p['P']-m['P'])/(2*d)),1e-9)

    def test_12_nonabelian_total_frame_invariance(self):
        q=mixed_static_EB(self.H,self.xy,self.P,frame='orbital',K=self.K)
        self.assertLess(abs(q['mixed_d2F_dEdB']-self.ref['mixed_d2F_dEdB']),1e-12)

    def test_13_nonabelian_omit_contact_fails(self):
        q=mixed_static_EB(self.H,self.xy,self.P,frame='orbital',K=self.K)
        self.assertGreater(abs(q['position_contact_dP_dB']),.01)
        self.assertGreater(abs(q['fermi_density_derivative']-self.ref['mixed_d2F_dEdB']),.01)

    def test_14_grand_fermi_E_hellmann_feynman(self):
        r=equilibrium_derivatives(self.H,self.xy,self.P);h=1e-5
        ep=equilibrium_derivatives(self.H,self.xy,self.P,E=h)['F']
        em=equilibrium_derivatives(self.H,self.xy,self.P,E=-h)['F']
        self.assertLess(abs((ep-em)/(2*h)-r['dF_dE']),1e-9)

    def test_15_grand_fermi_B_hellmann_feynman(self):
        r=equilibrium_derivatives(self.H,self.xy,self.P);h=1e-5
        p=equilibrium_derivatives(self.H,self.xy,self.P,replace(self.s,field=self.s.field+h))['F']
        m=equilibrium_derivatives(self.H,self.xy,self.P,replace(self.s,field=self.s.field-h))['F']
        self.assertLess(abs((p-m)/(2*h)-r['dF_dB']),1e-9)

    def test_16_grand_B_of_E_maxwell(self):
        x=finite_difference_mixed(self.H,self.xy,self.P,route='magnetic_of_electric')
        self.assertLess(abs(x-self.ref['mixed_d2F_dEdB']),2e-8)

    def test_17_grand_E_of_B_maxwell(self):
        x=finite_difference_mixed(self.H,self.xy,self.P,route='electric_of_magnetic')
        self.assertLess(abs(x-self.ref['mixed_d2F_dEdB']),2e-8)

    def test_18_grand_four_corner_energy(self):
        x=finite_difference_mixed(self.H,self.xy,self.P,route='free_energy')
        self.assertLess(abs(x-self.ref['mixed_d2F_dEdB']),1e-7)

    def test_19_landau_independent_Maxwell(self):
        x=finite_difference_mixed(self.H,self.xy,self.P,frame='landau_y',route='electric_of_magnetic')
        self.assertLess(abs(x-self.ref['mixed_d2F_dEdB']),2e-8)

    def test_20_nonabelian_independent_Maxwell(self):
        x=finite_difference_mixed(self.H,self.xy,self.P,frame='orbital',K=self.K)
        self.assertLess(abs(x-self.ref['mixed_d2F_dEdB']),2e-8)

    def test_21_fixed_mu_reference_value(self):
        self.assertAlmostEqual(self.ref['mixed_d2F_dEdB'],.028026623919872623,places=12)

    def test_22_fixed_N_reference_value(self):
        self.assertAlmostEqual(self.canon['mixed_d2F_dEdB_fixed_N'],.03877719176006167,places=12)

    def test_23_fixed_N_charge_Ward(self):
        self.assertLess(self.canon['particle_number_derivative_residual'],1e-12)
        self.assertGreater(abs(self.canon['mu_B_fixed_N']),.01)

    def test_24_fixed_N_B_of_E_independent_refit(self):
        x=finite_difference_canonical(self.H,self.xy,self.P,route='magnetic_of_electric')
        self.assertLess(abs(x-self.canon['mixed_d2F_dEdB_fixed_N']),1e-7)

    def test_25_fixed_N_E_of_B_independent_refit(self):
        x=finite_difference_canonical(self.H,self.xy,self.P,route='electric_of_magnetic')
        self.assertLess(abs(x-self.canon['mixed_d2F_dEdB_fixed_N']),1e-7)

    def test_26_fixed_N_Helmholtz_four_corner(self):
        x=finite_difference_canonical(self.H,self.xy,self.P,route='free_energy')
        self.assertLess(abs(x-self.canon['mixed_d2F_dEdB_fixed_N']),1e-7)

    def test_27_fixed_N_landau_invariance(self):
        x=canonical_mixed_static_EB(self.H,self.xy,self.P,frame='landau_y')
        self.assertLess(abs(x['mixed_d2F_dEdB_fixed_N']-self.canon['mixed_d2F_dEdB_fixed_N']),1e-12)

    def test_28_fixed_N_nonabelian_invariance(self):
        x=canonical_mixed_static_EB(self.H,self.xy,self.P,frame='orbital',K=self.K)
        self.assertLess(abs(x['mixed_d2F_dEdB_fixed_N']-self.canon['mixed_d2F_dEdB_fixed_N']),1e-12)

    def test_29_distinct_mu_vs_N(self):
        self.assertGreater(abs(self.ref['mixed_d2F_dEdB']-self.canon['mixed_d2F_dEdB_fixed_N']),.01)

    def test_30_constant_shift_fixed_N_invariant(self):
        x=canonical_mixed_static_EB(self.H,self.xy,self.P+2*np.eye(4))
        self.assertLess(abs(x['mixed_d2F_dEdB_fixed_N']-self.canon['mixed_d2F_dEdB_fixed_N']),1e-12)

    def test_31_constant_shift_fixed_mu_not_invariant(self):
        x=mixed_static_EB(self.H,self.xy,self.P+2*np.eye(4))
        self.assertGreater(abs(x['mixed_d2F_dEdB']-self.ref['mixed_d2F_dEdB']),.045)

    def test_32_real_TR_zero_B(self):
        h,xy,d=benchmark_graph(False)
        p=full_projected_polarization(xy,self.s.charge_over_hbar,dipoles=d)
        x=mixed_static_EB(h,xy,p,replace(self.s,field=0.))
        self.assertLess(abs(x['mixed_d2F_dEdB']),1e-13)

    def test_33_real_TR_B_oddness(self):
        h,xy,d=benchmark_graph(False)
        p=full_projected_polarization(xy,self.s.charge_over_hbar,dipoles=d)
        plus=mixed_static_EB(h,xy,p,replace(self.s,field=.23))['mixed_d2F_dEdB']
        minus=mixed_static_EB(h,xy,p,replace(self.s,field=-.23))['mixed_d2F_dEdB']
        self.assertLess(abs(plus+minus),1e-12)
        self.assertGreater(abs(plus),.002)

    def test_34_TR_broken_complex_loop_nonzero_B_zero(self):
        x=mixed_static_EB(self.H,self.xy,self.P,replace(self.s,field=0.))
        self.assertGreater(abs(x['mixed_d2F_dEdB']),.02)

    def test_35_degenerate_fermi_Frechet(self):
        h=np.diag([-.2,-.2,.9,1.4]).astype(complex)
        rng=np.random.default_rng(5194)
        x=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));m=(x+x.conj().T)/2
        d=1e-5
        def r(H):
            e,u=np.linalg.eigh(H)
            return (u*fermi(e,.1,.33))@u.conj().T
        finite=(r(h+d*m)-r(h-d*m))/(2*d)
        self.assertLess(np.linalg.norm(finite-finite_fermi_frechet(h,m,.1,.33)),1e-9)

    def test_36_reject_nonhermitian_position(self):
        with self.assertRaises(ValueError):
            full_projected_polarization(self.xy,-1,dipoles=np.eye(4)*1j)

    def test_37_reject_zero_T(self):
        with self.assertRaises(ValueError):Spec(temperature=0)

    def test_38_reject_bad_orbital_generator_and_frame(self):
        with self.assertRaises(ValueError):operators(self.H,self.xy,self.P,frame='orbital')
        with self.assertRaises(ValueError):operators(self.H,self.xy,self.P,frame='bad')

    def test_39_refuse_bad_finite_difference_steps(self):
        with self.assertRaises(ValueError):finite_difference_mixed(self.H,self.xy,self.P,step=0)
        with self.assertRaises(ValueError):finite_difference_canonical(self.H,self.xy,self.P,step=-1)

    def test_40_science_scope_not_transport(self):
        self.assertIn('NOT_DRIVEN',self.ref['status'])
        self.assertIn('NOT_TRANSPORT',self.canon['status'])


if __name__=='__main__':unittest.main()
