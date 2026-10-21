"""Original source SM Eq4-10 -> Eq20 adiabatic electric pole; no sign fit."""
import unittest
import numpy as np
from qgmap.sm_literal_spin_phe_v49 import SourceConfig,model,SX,SY,ID
from qgmap.sm_two_color_liouville_v50 import (
    Tone,_matrices,_resolvent,_rho0_source,two_color,run_limit)

class V50SMTwoColor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=SourceConfig(nr=80,nphi=104)
        cls.ref=model(cls.c,'main')['xy']['SM_Eq20_current_physical_B']

    def test_two_colours_are_independent_and_pdf_taus_not_global_eta(self):
        t=Tone(omega_a=.03,omega_b=-.02,eta_a=.005,eta_b=.012)
        self.assertNotEqual(t.eta_a,t.eta_b)
        a=two_color(SourceConfig(nr=24,nphi=40),t)
        self.assertEqual(a['tone']['omega_a'],.03)
        self.assertEqual(a['tone']['eta_b'],.012)
        self.assertIn('no all-sector RTA',a['source'])

    def test_band_resolvent_is_independent_matrix_liouville_solution(self):
        c=self.c;x=np.array([.27,-.31]);y=np.array([.36,.48])
        h,e,u,ud,rho=_matrices(x,y,c,'main')
        source=np.array([[[.12+.02j,.14+.13j],[-.22+.18j,.08+.02j]],
                         [[.01-.23j,.32+.19j],[.11+.12j,-.15+.07j]]])
        z=.003-1j*.052
        ans=_resolvent(h,e,u,ud,source,z)
        lhs=z*ans+1j*(h@ans-ans@h)
        np.testing.assert_allclose(lhs,source,atol=1e-13)

    def test_first_order_sm_eq6_and_eq7_in_original_eigenbasis(self):
        c=self.c;x=np.array([.27,-.31]);y=np.array([.36,.48])
        t=Tone(.031,-.043,.005,.009,1e-6)
        h,e,u,ud,rho=_matrices(x,y,c,'main')
        D0=_rho0_source(x,y,c,'main',t.dk,'y')
        s=ud@SX@u
        v=ud@(-SX)@u
        f=1/(1+np.exp(np.clip((e-c.mu)/c.temperature,-500,500)))
        za=t.eta_a-1j*t.omega_a;zb=t.eta_b-1j*t.omega_b
        e1=ud@_resolvent(h,e,u,ud,D0,za)@u
        b1=ud@_resolvent(h,e,u,ud,-1j*(SX@rho-rho@SX),zb)@u
        for j in range(len(x)):
            for m,n in ((0,1),(1,0)):
                eps=e[j,m]-e[j,n]
                r_mn=1j*v[j,m,n]/(e[j,n]-e[j,m])
                paper_E=(f[j,n]-f[j,m])*r_mn/(t.omega_a-eps+1j*t.eta_a)
                paper_B=(f[j,n]-f[j,m])*s[j,m,n]/(t.omega_b-eps+1j*t.eta_b)
                self.assertAlmostEqual(e1[j,m,n],paper_E,delta=5e-7)
                self.assertAlmostEqual(b1[j,m,n],paper_B,delta=1e-12)

    def test_full_two_frequency_dc_residue_matches_original_eq20_not_eq31(self):
        t=Tone(eta_a=.001,eta_b=.001)
        z=two_color(self.c,t,'main','xy')
        eq31=model(self.c,'main')['xy']['SM_printed_Eq31_current_physical_B']
        self.assertLess(abs(z['eta_a_times_total'].real-self.ref),3e-7)
        self.assertGreater(abs(z['eta_a_times_total'].real-eq31),.0049)
        self.assertLess(abs(z['eta_a_times_total'].imag),1e-8)
        self.assertEqual(z['scientific_status'],'PARTIAL')

    def test_two_individual_EB_orderings_large_cancel_only_after_sum(self):
        z=two_color(self.c,Tone(eta_a=.002,eta_b=.002))
        self.assertGreater(z['eta_a_times_b_then_e'].real,.003)
        self.assertLess(z['eta_a_times_e_then_b'].real,-.006)
        self.assertLess(z['eta_a_times_total'].real,0)
        self.assertLess(abs(z['eta_a_times_total'].real-self.ref),6e-7)
        self.assertGreater(abs(z['eta_a_times_total'].real-z['eta_a_times_e_then_b'].real),.003)

    def test_regulator_path_changes_individual_pieces_not_combined_pole(self):
        t1=two_color(self.c,Tone(eta_a=.002,eta_b=.0002))
        t2=two_color(self.c,Tone(eta_a=.0002,eta_b=.002))
        self.assertGreater(abs(t1['eta_a_times_b_then_e'].real-t2['eta_a_times_b_then_e'].real),.004)
        for z in (t1,t2):self.assertLess(abs(z['eta_a_times_total'].real-self.ref),4e-7)

    def test_both_main_and_SM_helicities_and_both_tensor_orderings(self):
        c=SourceConfig(nr=68,nphi=96)
        t=Tone(eta_a=.001,eta_b=.001,dk=2e-5)
        for direction in ('xy','yx'):
            a=two_color(c,t,'main',direction);b=two_color(c,t,'sm',direction)
            literal=model(c,'main')[direction]['SM_Eq20_current_physical_B']
            self.assertLess(abs(a['eta_a_times_total'].real-literal),1.2e-6)
            self.assertLess(abs(a['eta_a_times_total'].real+b['eta_a_times_total'].real),1e-10)

    def test_eta_a_pole_is_separate_from_finite_omega_a(self):
        c=SourceConfig(nr=48,nphi=64)
        a=two_color(c,Tone(omega_a=.04,eta_a=.002,eta_b=.001))
        b=two_color(c,Tone(omega_a=.04,eta_a=.001,eta_b=.001))
        self.assertLess(abs(a['eta_a_times_total']),.0003)
        self.assertLess(abs(b['eta_a_times_total']),abs(a['eta_a_times_total'])*.55)
        with self.assertRaises(ValueError):two_color(c,Tone(omega_a=.04),pole_only=True)

    def test_finite_frequency_complex_conjugate_real_time_reciprocity(self):
        c=SourceConfig(nr=36,nphi=56)
        for wa,wb in ((0.,.07),(.03,-.07)):
            a=two_color(c,Tone(omega_a=wa,omega_b=wb,eta_a=.004,eta_b=.002))
            b=two_color(c,Tone(omega_a=-wa,omega_b=-wb,eta_a=.004,eta_b=.002))
            self.assertAlmostEqual(a['physical_B_total'],b['physical_B_total'].conjugate(),delta=2e-9)
            self.assertGreater(abs(a['physical_B_total'].imag),1e-4)

    def test_explicit_electric_population_residue_matches_SM_eq20(self):
        t=Tone(eta_a=.001,eta_b=.001)
        a=two_color(self.c,t,pole_only=True)
        v=a['eta_a_times_physical_B_current_residue']
        self.assertLess(abs(v.real-self.ref),3e-7)
        self.assertLess(abs(v.imag),1e-11)

    def test_independent_k_derivative_precision(self):
        c=SourceConfig(nr=48,nphi=64)
        p=two_color(c,Tone(eta_a=.003,eta_b=.002,dk=2e-5))
        q=two_color(c,Tone(eta_a=.003,eta_b=.002,dk=6e-5))
        self.assertLess(abs(p['eta_a_times_total']-q['eta_a_times_total']),2e-7)

    def test_fail_closed_and_original_source_sign_preserved(self):
        with self.assertRaises(ValueError):Tone(eta_a=0)
        with self.assertRaises(ValueError):Tone(eta_b=-.1)
        with self.assertRaises(ValueError):two_color(self.c,helicity='other')
        with self.assertRaises(ValueError):two_color(self.c,direction='zz')
        with self.assertRaises(ValueError):run_limit(self.c,eta_as=())
        self.assertLess(self.ref,0.)
        self.assertGreater(model(self.c,'main')['xy']['SM_printed_Eq31_current_physical_B'],0)

if __name__=='__main__':unittest.main()
