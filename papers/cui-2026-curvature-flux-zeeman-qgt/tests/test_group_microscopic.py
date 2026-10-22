"""Round16: 31 microscopic four-band MPG models, 62 Jahn ranks, causal guards."""
import unittest
import numpy as np
from flux_lab.group_microscopic import (
    build_group,momentum_map,spin_axial,symmetrized_model,
    local_group_residual,sector_nullspace,verify_catalog_projectors,
    project_frequency_tensor,clean_bz_even_odd)
from flux_lab.magnetic_symmetry import load_groups
from flux_lab.causal_response import interband_kubo
K=np.array([.41,-.32]); Z=.051+.011j
ALL=load_groups()

class GroupMicroscopicBridge(unittest.TestCase):
    def test_31_magnetic_point_group_spinors_close(self):
        self.assertEqual(len(ALL),31)
        for name,g in ALL.items():
            self.assertEqual(len(build_group(name)),len(g["operations"]))

    def test_actual_tr_momentum_is_negative_not_realspace_identity(self):
        g=next(q for q in build_group("11'") if q.anti)
        np.testing.assert_allclose(g.R,np.eye(2),atol=1e-13)
        np.testing.assert_allclose(momentum_map(g),-np.eye(2),atol=1e-13)
        np.testing.assert_allclose(spin_axial(g),-np.eye(3),atol=1e-13)

    def test_c2zt_momentum_fixed_but_realspace_rotates(self):
        g=next(q for q in build_group("2'") if q.anti)
        np.testing.assert_allclose(g.R,-np.eye(2),atol=1e-13)
        np.testing.assert_allclose(momentum_map(g),np.eye(2),atol=1e-13)
        np.testing.assert_allclose(spin_axial(g),np.diag([1,1,-1]),atol=1e-13)

    def test_mirror_tr_momentum_and_spin_signs(self):
        g=next(q for q in build_group("m'") if q.anti)
        np.testing.assert_allclose(g.R,np.diag([-1,1]),atol=1e-13)
        np.testing.assert_allclose(momentum_map(g),np.diag([1,-1]),atol=1e-13)
        np.testing.assert_allclose(spin_axial(g),np.diag([-1,1,1]),atol=1e-13)

    def test_microscopic_rank_equals_independent_jahn_for_all_31(self):
        rows=verify_catalog_projectors()
        self.assertEqual(len(rows),62)
        self.assertTrue(all(item["agree"] for item in rows),rows)

    def test_distinct_gray_and_rotation_antigroup_ranks(self):
        self.assertEqual([sector_nullspace("11'",p) for p in ("even","odd")],[4,0])
        self.assertEqual([sector_nullspace("4'",p) for p in ("even","odd")],[2,2])
        self.assertEqual([sector_nullspace("m'",p) for p in ("even","odd")],[2,2])

    def test_all_31_nondegenerate_fourband_kubo_covariance(self):
        for name in ALL:
            with self.subTest(group=name):
                r=local_group_residual(name,K,z=Z)
                self.assertLess(max(r[x] for x in ("H","velocity","chi","alpha")),1e-12,r)

    def test_independent_orbital_liouville_across_seven_groups(self):
        for name in ("11'","2'","m'","m'm'2","4'","4mm1'","6m'm'"):
            with self.subTest(group=name):
                r=local_group_residual(name,K,z=Z,liouville=True)
                self.assertLess(max(r[x] for x in ("H","velocity","chi","alpha")),2e-11,r)

    def test_analytic_velocity_matches_finite_difference(self):
        for name in ("11'","m'","4'","4mm"):
            h,v,s=symmetrized_model(K,name)
            for a in range(2):
                dx=np.eye(2)[a]*2e-6
                hp=symmetrized_model(K+dx,name)[0]
                hm=symmetrized_model(K-dx,name)[0]
                np.testing.assert_allclose(v[a],(hp-hm)/(4e-6),atol=1e-9)

    def test_group_frequency_projection_is_idempotent(self):
        X=np.array([[.23+.1j,-.3,.4],[.57,.11-.2j,-.8]])
        for name in ("11'","m'","4'","6'","6mm"):
            for parity in ("even","odd"):
                y=project_frequency_tensor(name,X,parity)
                np.testing.assert_allclose(y,project_frequency_tensor(name,y,parity),atol=5e-13)

    def test_square_bz_exact_symmetry_channels(self):
        for name in ("11'","2'","m'","4'","4mm"):
            even,odd=clean_bz_even_odd(name,n=8)
            np.testing.assert_allclose(even,project_frequency_tensor(name,even,"even"),atol=3e-12)
            np.testing.assert_allclose(odd,project_frequency_tensor(name,odd,"odd"),atol=3e-12)

    def test_tr_odd_cancels_even_survives(self):
        even,odd=clean_bz_even_odd("11'",n=8)
        self.assertLess(np.max(abs(odd)),1e-13)
        self.assertGreater(np.max(abs(even)),1e-7)

    def test_allowed_rank_need_not_produce_nonzero_sample(self):
        self.assertEqual(sector_nullspace("2'","even"),4)
        even,odd=clean_bz_even_odd("2'",n=8)
        self.assertLess(np.max(abs(even[:,:2])),1e-12)
        self.assertLess(np.max(abs(odd[:,:2])),1e-12)

    def test_anti_mirror_allowed_odd_response_nonzero(self):
        even,odd=clean_bz_even_odd("m'",n=8)
        self.assertGreater(abs(odd[1,1]),5e-5)

    def test_wrong_realspace_momentum_mapping_negative_control(self):
        g=next(x for x in build_group("11'") if x.anti)
        H,V,S=symmetrized_model(K,"11'")
        U4=np.kron(np.eye(2),g.U)
        t=lambda A:U4@A.conj()@U4.conj().T
        self.assertGreater(np.max(abs(H-t(H))),.1)
        np.testing.assert_allclose(symmetrized_model(momentum_map(g)@K,"11'")[0],t(H),atol=1e-12)

    def test_wrong_retarded_instead_of_advanced_negative_control(self):
        for name in ("2'","m'"):
            g=next(x for x in build_group(name) if x.anti)
            H,V,S=symmetrized_model(K,name)
            HP,VP,SP=symmetrized_model(momentum_map(g)@K,name)
            z=.051+.11j
            chi=interband_kubo(HP,VP,SP,z=z,mu=.17,temperature=.17)[0]
            wrong=interband_kubo(H,V,S,z=-z.real+1j*z.imag,mu=.17,temperature=.17)[0]
            mismatch=np.max(abs(chi-momentum_map(g)@wrong@spin_axial(g).T))
            self.assertGreater(mismatch,1e-5)

    def test_hexagonal_reciprocal_bz_not_square_grid(self):
        with self.assertRaisesRegex(ValueError,"hexagonal"):
            clean_bz_even_odd("6",n=8)

    def test_invalid_group_and_sector_inputs_fail_closed(self):
        with self.assertRaises(ValueError):sector_nullspace("1","invalid")
        with self.assertRaises(ValueError):project_frequency_tensor("1",np.eye(2),"even")
        with self.assertRaises(ValueError):build_group("unknown")
if __name__=="__main__":unittest.main()
