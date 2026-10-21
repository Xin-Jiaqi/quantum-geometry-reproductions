"""v51: source-index algebra, transformation law, spin Ward Eq(33) independent gates."""
import unittest
import numpy as np
from qgmap.sm_sign_covariance_v51 import (
    I,SX,SY,SZ,source_sign_identity,hamiltonian,single_point,
    finite_difference_spin_ward,integrated_table)
from qgmap.sm_literal_spin_phe_v49 import SourceConfig,geometry
from qgmap.sm_two_color_liouville_v50 import Tone,two_color


class SMSpinSignCovarianceV51(unittest.TestCase):
    def test_arbitrary_hermitian_two_band_matrix_algebra_not_specific_dirac(self):
        rng=np.random.default_rng(1151)
        for _ in range(120):
            v=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            s=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            vv=(v+v.conj().T)/2;ss=(s+s.conj().T)/2
            df=rng.normal(size=2)
            z=source_sign_identity(np.array([-.31,.78]),vv,ss,df)
            self.assertLess(abs(z['direct_minus_Z']),2e-12)
            self.assertLess(abs(z['eq20_plus_printed_eq31']),2e-12)
            self.assertLess(abs(z['imaginary_spurious']),2e-12)

    def test_independent_band_U1_gauge_leaves_two_branches_and_relative_sign(self):
        r=np.random.default_rng(42)
        for i in range(20):
            vs=r.normal(size=(2,2))+1j*r.normal(size=(2,2))
            sp=r.normal(size=(2,2))+1j*r.normal(size=(2,2))
            vs=(vs+vs.conj().T)/2;sp=(sp+sp.conj().T)/2
            df=np.array([-.16,.35]);e=np.array([-.9,.28])
            original=source_sign_identity(e,vs,sp,df)
            phases=np.exp(1j*r.uniform(-np.pi,np.pi,size=2))
            d=np.diag(phases)
            new=source_sign_identity(e,d.conj().T@vs@d,d.conj().T@sp@d,df)
            self.assertAlmostEqual(original['j20_per_Bbar'],new['j20_per_Bbar'],delta=1e-12)
            self.assertAlmostEqual(original['printed_eq31_per_Bbar'],new['printed_eq31_per_Bbar'],delta=1e-12)

    def test_mixed_operator_antiunitary_not_needed_spin_SU2_is_enough(self):
        for x,y in ((.4,.2),(-.31,.62)):
            h=hamiltonian(x,y,helicity='main')
            np.testing.assert_allclose(SZ@h@SZ,hamiltonian(x,y,helicity='sm'),atol=1e-14)
            np.testing.assert_allclose(SZ@SX@SZ,-SX,atol=1e-14)
            np.testing.assert_allclose(SZ@SY@SZ,-SY,atol=1e-14)
            np.testing.assert_allclose(SZ@(0.125*I+SY)@SZ,0.125*I-SY,atol=1e-14)

    def test_representations_transform_same_physical_spin_response(self):
        t=integrated_table()
        self.assertEqual(t['scientific_status'],'PARTIAL')
        for v in t['basis_covariance_residuals'].values():
            for err in v.values():self.assertLess(abs(err),1e-13)
        self.assertGreater(t['table']['main_spin_x']['sym_eq20'],0)
        self.assertLess(t['table']['sm_spin_x']['sym_eq20'],0)
        self.assertAlmostEqual(t['table']['main_spin_x']['sym_eq20'],
                               -t['table']['sm_spin_x']['sym_eq20'],delta=1e-12)

    def test_eq20_31_opposite_sign_in_all_four_spin_representation_choices(self):
        t=integrated_table()
        for d in t['table'].values():
            self.assertAlmostEqual(d['sym_eq20'],-d['sym_printed_eq31'],delta=1e-12)
            for ab in ('xy','yx'):
                self.assertAlmostEqual(d[ab]['eq20_per_Bbar'],
                                       -d[ab]['printed_eq31_per_Bbar'],delta=1e-12)

    def test_eq33_local_completeness_is_separately_valid(self):
        for helicity in ('main','sm'):
            for ssign in (1,-1):
                result=finite_difference_spin_ward(helicity=helicity,spin_sign=ssign)
                for alpha in ('x','y'):
                    self.assertLess(result[alpha]['maximum_residual'],2e-9)

    def test_dirac_local_Z_is_derivative_of_spin_texture(self):
        x,y=.29,.41;rr=np.hypot(x,y)
        values=finite_difference_spin_ward(x,y)
        self.assertAlmostEqual(values['x']['Z_band'][1],x*y/rr**3,delta=1e-12)
        self.assertAlmostEqual(values['y']['Z_band'][1],-x*x/rr**3,delta=1e-12)

    def test_eq33_integrated_BZ_ibp_symmetry_is_independent_of_eq31_sign(self):
        t=integrated_table(nr=120,nphi=160)
        for d in t['table'].values():
            self.assertLess(abs(d['xy']['eq20_per_Bbar']-d['yx']['eq20_per_Bbar']),2e-8)
        self.assertGreater(abs(t['table']['main_spin_x']['sym_eq20']),.002)

    def test_spin_basis_flip_does_not_flip_relative_ratio_for_beta_sweep(self):
        for beta in (.025,.125,.3,.6):
            t=integrated_table(beta=beta,nr=84,nphi=120)
            for d in t['table'].values():
                self.assertAlmostEqual(d['sym_eq20'],-d['sym_printed_eq31'],delta=3e-12)

    def test_mu_sign_and_spin_sign_are_different_transformations(self):
        p=integrated_table(mu=.3,nr=90,nphi=120)
        n=integrated_table(mu=-.3,nr=90,nphi=120)
        for name in p['table']:
            self.assertAlmostEqual(p['table'][name]['sym_eq20'],
                                   -n['table'][name]['sym_eq20'],delta=5e-7)

    def test_v50_source_liouville_sign_is_same_as_eq20_under_field_basis_transform(self):
        c=SourceConfig(nr=44,nphi=64)
        tone=Tone(eta_a=.001,eta_b=.001)
        x=two_color(c,tone,helicity='main',direction='xy')
        y=two_color(c,tone,helicity='sm',direction='xy')
        self.assertLess(x['eta_a_times_total'].real,0)
        self.assertGreater(y['eta_a_times_total'].real,0)
        self.assertAlmostEqual(x['eta_a_times_total'].real,-y['eta_a_times_total'].real,delta=1e-9)

    def test_reject_nonhermitian_or_degenerate_input(self):
        e=np.array([0.,1.]);v=np.array([[0,1],[1,0]],complex)
        with self.assertRaises(ValueError):source_sign_identity([0,0],v,SX,[-.2,.3])
        with self.assertRaises(ValueError):source_sign_identity(e,np.array([[0,1],[0,0]],complex),SX,[-.2,.3])
        with self.assertRaises(ValueError):source_sign_identity(e,v,SX,[-.2])
        with self.assertRaises(ValueError):hamiltonian(.2,.4,helicity='wrong')
        with self.assertRaises(ValueError):single_point(0.,0.)
        with self.assertRaises(ValueError):finite_difference_spin_ward(step=-1)


if __name__=='__main__':unittest.main()
